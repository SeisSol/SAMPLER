#=
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-07-13
=#

using Interpolations
using Base.Threads
using Printf
using FFTW
using NetCDF
using Base.Filesystem

const G = begin
    #              ∞     (-1)^n * (2n + 1)
    # G(r) ≈ 1/π * Σ ——————————————————————————
    #             n=0 [(2n + 1)^2 + r^2]^(3/2)
    #
    # r ∈ [0; √(800)]
    Δr = .01
    l_r = 0:Δr:(sqrt(800) + Δr)

    function Σ(r::Float64)
        Σ = 0.
        for n ∈ 10000:-1:0
            frac = (2n + 1) / ((2n + 1)^2 + r^2)^1.5
            if n & 0x1 != 0 # n is odd, subtract
                Σ -= frac
            else # n is even, add
                Σ += frac
            end
        end

        return Σ
    end

    l_G = [1/Float64(π) * Σ(r) for r in l_r]

    itp = interpolate(l_G, BSpline(Constant()))
    extrapolate(scale(itp, l_r), Flat())
end

function precalculate_σ(h_min :: Float64, h_max :: Float64, Δx :: Float64, Δy :: Float64; n_h :: Float64 = 20.)
    Δh = (h_max - h_min) / 10
    if Δh == 0.; Δh = 1.; end
    l_h = max(0.0001, h_min):Δh:(h_max + Δh)

    #                       h^2
    # σ = ———————————————————————————————————————
    #       Σ   Σ  G(√[(nΔx)^2 + (mΔy)^2]/h)ΔxΔy
    #      n∈N m∈M
    #
    # N = {n ∈ Z| |nΔx| ≤ n_h * h}
    # M = {m ∈ Z| |mΔy| ≤ n_h * h}
    #
    # As the domain of n and m are point symmetric around 0, only ~1/4th of the sum element have to be evaluated:
    # * One quadrant (n > 0, m > 0), multiplied by 4
    # * One half of the x-axis (n > 0), multiplied by 2
    # * One hald of the y-axis (m > 0), multiplied by 2
    # * The point at the origin (n, m) = (0, 0)

    function σ(h)
        σ_inv = 0.
        n_max = ceil(Int, (n_h * h / Δx))
        m_max = ceil(Int, (n_h * h / Δy))

        N = n_max:-1:1
        M = m_max:-1:1

        for n ∈ N, m ∈ M
            σ_inv += 4*G(sqrt((n*Δx)^2 + (m*Δy)^2) / h)
        end

        for n ∈ N
            σ_inv += 2*G(n*Δx / h)
        end

        for m ∈ M
            σ_inv += 2*G(m*Δy / h)
        end

        σ_inv += G(0)

        return h^2 / (σ_inv * Δx * Δy)
    end

    l_σ = [σ(h) for h ∈ l_h]

    itp = interpolate(l_σ, BSpline(Cubic(Flat(OnGrid()))))
    return extrapolate(scale(itp, l_h), Flat())
end

function apply_kajiura!(b::AbstractArray{Float64, 2}, d::AbstractArray{Float64, 2}, η::AbstractArray{Float64, 2}, 
                h_min :: Float64, h_max :: Float64, Δx :: Float64, Δy :: Float64; water_level :: Float64 = 0., n_h :: Float64 = 20., 
                σ = precalculate_σ(h_min, h_max, Δx, Δy; n_h=n_h))
    
    filter_nx_half = ceil(Int, n_h * h_max / Δx / 2)
    filter_ny_half = ceil(Int, n_h * h_max / Δy / 2)

    nx = size(η, 2)
    ny = size(η, 1)

    filter = zeros(Float64, (ny, nx))

    @assert filter_nx_half < nx
    @assert filter_ny_half < ny

    println("  Computing filter matrix...")
    Threads.@threads for x ∈ -filter_nx_half:filter_nx_half
        for y ∈ -filter_ny_half:filter_ny_half
            filter[((ny + y) % ny) + 1, ((nx + x) % nx) + 1] = G(sqrt((x * Δx)^2 + (y * Δy)^2) / h_max)
        end
    end

    println("  Computing displacement matrix...")
    Threads.@threads for x ∈ 1:size(η, 2)
        for y ∈ 1:size(η, 1)
            h_yx = max(0., water_level - b[y, x]) # height 0 on land

            η[y, x] = 
                if h_yx != 0.
                    σ(h_yx) * Δx * Δy / h_yx^2 * d[y, x]
                else
                    0. # No displacement where no water is (Kajiura not applicable)
                end
        end
    end

    println("  Computing filter FFT...")
    Filter = fft(filter)

    println("  Computing displacement FFT...")
    Η = fft(η)

    println("  Computing product...")
    Η .*= Filter

    println("  Computing IFFT...")
    η_complex = reinterpret(Float64, ifft(Η))
    copyto!(η, @view η_complex[1:2:end-1])
end

function main()
    if length(ARGS) != 3
        println("Usage: julia ./kajiura.jl <in_file.nc> <out_file.nc> <timestep_end>")
        return
    end

    println("Using $(nthreads()) threads.")

    in_filename = ARGS[1]
    out_filename = ARGS[2]
    t_end = parse(Int, ARGS[3])
    nc = NetCDF.open(in_filename, mode=NC_NOWRITE)
    d = nc["d"]

    lx = nc["x"]
    ly = nc["y"]

    Δx = lx[2] - lx[1]
    Δy = ly[2] - ly[1]

    nx = length(lx)
    ny = length(ly)

    l_times = nc["time"]
    n_times = length(l_times)
    last_timestamp = l_times[end]

    if t_end == -1 || t_end > n_times - 1
	    t_end = n_times - 1
    end

    println("Processing $nx × $ny cells over $(t_end + 1) timesteps.")

    b = Array{Float64, 2}(undef, (ny, nx))
    ncread!(in_filename, "b", b)
    current_disp = zeros(ny, nx)
    current_η_diff = Array{Float64, 2}(undef, (ny, nx))
    current_d_diff = Array{Float64, 2}(undef, (ny, nx))

    println("Creating output file")

    isfile(out_filename) && rm(out_filename)    
    timeatts = Dict("units" => "seconds")
    xatts = Dict("units" => "m")
    yatts = Dict("units" => "m")

    nccreate(out_filename, "eta_diff", "y", ly, yatts, "x", lx, xatts, "time", l_times[1:t_end+1], timeatts)
    nccreate(out_filename, "d_diff", "y", "x", "time")

    nccreate(out_filename, "b", "y", "x")
    ncwrite(b, out_filename, "b")

    for t in 1:t_end + 1
        println("Working on timestep $(t - 1) of $(t_end)")
        current_η_diff .= 0.
        current_d_diff .= d[:,:,t] .- current_disp
        apply_kajiura!(b .+ current_disp, current_d_diff, current_η_diff, -maximum(b), -minimum(b), Δx, Δy)
        current_disp = d[:,:,t]

        println("  Writing output for timestep")
        ncwrite(current_d_diff, out_filename, "d_diff", start=[1,1,t], count=[-1,-1,1])
        ncwrite(current_η_diff, out_filename, "eta_diff", start=[1,1,t], count=[-1,-1,1])
    end

    println("Done.")
end

main()

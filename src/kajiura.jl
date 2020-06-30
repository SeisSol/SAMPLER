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
    # r ∈ [0; √(50)]
    Δr = .01
    l_r = 0:Δr:(sqrt(50) + Δr)

    function Σ(r::Float64)
        Σ = 0.
        for n ∈ 100:-1:0
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

function precalculate_σ(h_min :: Float64, h_max :: Float64, Δx :: Float64, Δy :: Float64; n_h :: Float64 = 5.)
    Δh = (h_max - h_min) / 10
    if Δh == 0.; Δh = 1.; end
    l_h = h_min:Δh:(h_max + Δh)

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
                h_min :: Float64, h_max :: Float64, Δx :: Float64, Δy :: Float64; water_level :: Float64 = 0., n_h :: Float64 = 5., 
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
            h_yx = water_level - b[y, x]
            η[y, x] = σ(h_yx) * Δx * Δy / h_yx^2 * d[y, x]
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

    println("  Done!")
end

function main()
    if length(ARGS) != 2
        println("Usage: julia ./kajiura.jl <in_file.nc> <out_file.nc>")
        exit(0)
    end

    println("Using $(nthreads()) threads.")

    in_filename = ARGS[1]
    out_filename = ARGS[2]
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

    println("Processing $nx × $ny cells over $n_times timesteps.")

    b = Array{Float64, 2}(undef, (ny, nx))
    ncread!(in_filename, "b", b)
    current_disp = zeros(ny, nx)
    current_η = zeros(ny, nx)

    println("Creating output file")

    isfile(out_filename) && rm(out_filename)    
    nccreate(out_filename, "h", "y", ly, "x", lx)
    buffer_η = Array{Float64, 2}(undef, (ny, nx))

    for t in 1:n_times
        println("Working on timestep $t of $n_times")
        buffer_η .* 0.
        apply_kajiura!(b .+ current_disp, d[:,:,t] .- current_disp, buffer_η, -minimum(b), -maximum(b), Δx, Δy)
        current_disp = d[:,:,t]
        current_η .+= buffer_η
    end

    println("Writing outputs...")

    nccreate(out_filename, "b", "y", "x")
    nccreate(out_filename, "d", "y", "x")

    ncwrite(b .+ current_disp, out_filename, "b")
    ncwrite(zeros(Float64, (ny, nx)), out_filename, "d")
    ncwrite(current_η, out_filename, "h")

    println("Done.")
end

main()
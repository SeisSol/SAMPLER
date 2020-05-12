module Kajiura
    using Interpolations
    using Base.Threads
    using Printf

    export apply_kajiura!, precalculate_σ

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
        l_h = h_min:Δh:(h_max + Δh)

        #                     h^2
        # σ = ———————————————————————————————————
        #       Σ   Σ  G(√[(nΔx)^2 + (mΔy)^2]/h)
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
                   h_min :: Float64, h_max :: Float64, Δx :: Float64, Δy :: Float64, water_level :: Float64; n_h :: Float64 = 5., 
                   σ = precalculate_σ(h_min, h_max, Δx, Δy; n_h=n_h))
        
        ny, nx = size(b)

        n_max = ceil(Int, n_h * h_max / Δx)
        m_max = ceil(Int, n_h * h_max / Δy)

        println("h_max=$h_max, h_min=$h_min, n_max=$n_max, m_max=$m_max")

        mmn_quantities = Array{NTuple{2, Float64}, 2}(undef, (ny, nx))
        mji_abs = Array{Float64, 2}(undef, (m_max+1, n_max+1))

        i_nonzero_x = [[nx, 1] for i ∈ 1:nthreads()]
        i_nonzero_y = [[ny, 1] for i ∈ 1:nthreads()]
        l_nonzero_cells = zeros(Int, nx)

        #============================================#
        # Precalculate factors that are constant for
        # a specific bottom displacement D_ij.
        # As these will be used in 
        #============================================#

        println("h_max_b=", water_level-minimum(b), "; h_min_b=", water_level-maximum(b))

        Threads.@threads for i ∈ 1:nx
            tid = threadid()

            for j ∈ 1:ny
                D_ij = d[j, i]
                if abs(D_ij) < eps()
                    mmn_quantities[j, i] = (0., NaN64)
                    continue
                end

                h_ij = water_level-b[j, i]
                σ_ij = σ(h_ij)
                mmn_quantities[j, i] = (Δx * Δy / h_ij^2 * σ_ij * D_ij, h_ij)
                l_nonzero_cells[i] += 1

                # TODO thread safety
                if i_nonzero_x[tid][1] > i; i_nonzero_x[tid][1] = i end
                if i_nonzero_y[tid][1] > j; i_nonzero_y[tid][1] = j end
                if i_nonzero_x[tid][2] < i; i_nonzero_x[tid][2] = i end
                if i_nonzero_y[tid][2] < j; i_nonzero_y[tid][2] = j end
            end
        end

        i_nonzero_x = (minimum(rng -> rng[1], i_nonzero_x), maximum(rng -> rng[2], i_nonzero_x))
        i_nonzero_y = (minimum(rng -> rng[1], i_nonzero_y), maximum(rng -> rng[2], i_nonzero_y))

        if i_nonzero_x[1] > i_nonzero_x[2]
            println("  No non-zero displacements found, skipping Kajiura filter.")
            return
        end

        Threads.@threads for i ∈ 0:n_max
            for j ∈ 0:m_max
                mji_abs[j+1, i+1] = sqrt((i * Δx)^2 + (j*Δy)^2)
            end
        end

        n_rng = max(1, i_nonzero_x[1]-n_max):min(nx, i_nonzero_x[2]+n_max)
        m_rng = max(1, i_nonzero_y[1]-m_max):min(ny, i_nonzero_y[2]+m_max)

        @printf("  Applying Kajiura filter to %.2f%% of the domain.\n", length(n_rng)*length(m_rng)/(nx*ny)*100)
        println("  Kajiura: Kernel size is up to ", (2*n_max+1) ,"×", (2*m_max+1))

        n_threads = nthreads()
        n_nonzero_cells = sum(l_nonzero_cells)
        n_cells_per_thread = round(Int, n_nonzero_cells / n_threads)
        l_start_cols = Array{Int, 1}(undef, n_threads)

        l_start_cols[1] = n_rng[1]
        for i ∈ 2:n_threads
            start_col = l_start_cols[i-1]
            n_cells = 0
            while n_cells < n_cells_per_thread && start_col != nx
                n_cells += l_nonzero_cells[start_col]
                start_col += 1
            end
            l_start_cols[i] = start_col
        end

        Threads.@threads for tid ∈ 1:n_threads
            n_end_col = (tid == n_threads) ? n_rng[end] : (l_start_cols[tid+1] - 1)
            thread_n_rng = l_start_cols[tid]:n_end_col
            if thread_n_rng[1] > nx || n_end_col > nx; continue end

            cols_processed = 0
            cols_total = length(thread_n_rng)

            for n ∈ thread_n_rng
                for m ∈ m_rng
                    for i ∈ max(1, n - n_max):min(nx, n + n_max), 
                        j ∈ max(1, m - m_max):min(ny, m + m_max)

                        @inbounds q = mmn_quantities[j, i] # = (factor_ij, h_ij)
                        if abs(q[1]) < eps(); continue end
                        @inbounds η[m, n] += q[1] * G(mji_abs[abs(m-j)+1, abs(n-i)+1]/q[2])
                    end
                end
                if threadid() == n_threads ÷ 2
                    cols_processed += 1
                    @printf("  Kajiura: %.2f%% done.\n", (cols_processed)/cols_total*100)
                end
            end
        end
    end
end
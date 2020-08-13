#=
rasterization.jl:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

module Rasterization
    using Base.Threads
    using LinearAlgebra
    using Dates
    using Printf
    using NetCDF
    using Profile
    using Main.Util
    using Main.NC

    export rasterize, ZRange

    @enum ZRange z_all=-1 z_floor=0 z_surface=1
    @enum LoadBalancer naive=0 count=1 workload=2

    struct HesseNormalForm
        n0  :: NTuple{3, Float64}
        d   :: Float64
    end

    const EDGE_TOL = .0      # Tolerance for whether to consider points that are close to a simplex edge/face to be on that edge/face
    const Z_RANGE_TOL = .001 # Tolerance for points close to the z_range-limit. Points over the limit by less than Z_RANGE_TOL are not discarded

    """
        Rasterize variables stored in an unstructured grid into a regular grid. Write the results to a NetCDF file.

    # Arguments
    - `simplices`:      Each column consists of 3 or 4 point indices representing a triangle/tetrahedron
    - `points`:         Each column consists of an x-, y- and z-coordinate representing a point in 3D space
    - `in_vars`:        A 2D array of arrays. Each column corresponds to one variable and contains `size(times)` arrays 
                        and each of those arrays contains the variable value for each simplex (in the same order)
    - `in_var_names`:   The names, of the variables in `in_vars` in the same order as those vars
    - `times`:          An array of timestamps associated to the timesteps
    - `sampling_rate`:  A tuple (Δx, Δy, Δz) defining the regular grid cell size.
    - `out_filename`:   The path/name of the output NetCDF file
    - `mem_limit`:      The soft memory (in bytes) the software is allowed to utilize (soft limit!)
    - `z_range`:        A filter for discarding simplices above/below a threshold. `z_floor` discards everything above `0-Z_RANGE_TOL`, 
                        `z_surface` discards everything below. `z_all` keeps all simplices.
    - `t_begin`:        The first timestep to be rasterized
    - `t_end`:          The last timestep to be rasterized
    - `create_file`:    If `true`: create output file (already existing file will be overwritten). Else: use existing output file.
    - `load_balancer`:  The load balancer to be used.
    - `lb_params`:      The parameters for the workload-LB as obtained by the autotuning procedure
    - `lb_autotune`:    Whether to run the autotune procedure for the workload-LB. Only runs the first iteration of the rasterization.
    - `water_height`:   The offset to translate the water level to be at `z=0`.
    """
    function rasterize(
        simplices       :: AbstractArray{INDEX_TYPE, 2}, 
        points          :: AbstractArray{Float64, 2},
        in_vars         :: AbstractArray{A where A <: AbstractArray, 2},
        in_var_names    :: AbstractArray,
        times           :: AbstractArray{Float64, 1}, 
        sampling_rate   :: NTuple{3, Float64}, 
        out_filename    :: AbstractString,
        mem_limit       :: Int64; 
        z_range         :: ZRange = z_all, 
        t_begin         :: Integer = 1, 
        t_end           :: Integer = length(times), 
        create_file     :: Bool = false, 
        load_balancer   :: LoadBalancer = naive,
        lb_params       :: NTuple{2, Float64} = (0., 0.),
        lb_autotune     :: Bool = false,
        water_height    :: Float64 = 0.) where INDEX_TYPE <: Integer

        #============================================#
        # Code style & conventions
        #============================================#
        #
        # Some variables are prefixed with <prefix>_
        # to indicate their type or unit:
        #   i_  : Interval / Pair / 2-Tuple
        #   n_  : Count / Number of Objects
        #   o_  : Array Offset
        #   b_  : Size in Bytes
        #   l_  : 1D-Array / List
        #   d_  : Date & Time
        #   t_  : Timestep Index
        #   mMN_: M×N Matrix
        #   vN_ : N-Dimensional Vector
        #   s_  : String
        #============================================#

        i_domain_x = (minimum(points[1,:]), maximum(points[1,:]))
        i_domain_y = (minimum(points[2,:]), maximum(points[2,:]))
        i_domain_z = (minimum(points[3,:]), maximum(points[3,:]))

        println("Domain is $i_domain_x × $i_domain_y × $i_domain_z.")

        #============================================#
        # Set z-range filter
        #============================================#

        if z_range != z_all
            s_z_range = nothing

            if z_range == z_floor
                i_domain_z = (i_domain_z[1] - Z_RANGE_TOL, i_domain_z[2] - Z_RANGE_TOL)
                s_region = "floor"
            else # if z == z_surface
                i_domain_z = (i_domain_z[2] - Z_RANGE_TOL, i_domain_z[2] + Z_RANGE_TOL)
                s_region = "surface"
            end

            println("Only processing triangles on the ", s_region, " (z-domain limited to ", i_domain_z, ").")
        end

        n_samples_x = ceil(Int, (i_domain_x[end] - i_domain_x[1]) / sampling_rate[1])
        n_samples_y = ceil(Int, (i_domain_y[end] - i_domain_y[1]) / sampling_rate[2])
        n_samples_z = ceil(Int, (i_domain_z[end] - i_domain_z[1]) / sampling_rate[3])

        println("Sampling into $n_samples_x × $n_samples_y × $n_samples_z cells.")

        #============================================#
        # Calculate problem size & hardware limitations
        #============================================#

        n_threads                       = nthreads()

        n_simplex_points                = size(simplices, 1)
        n_dims                          = n_simplex_points - 1
        n_simplices                     = size(simplices, 2)
        n_simplices_per_thread          = n_simplices / n_threads

        b_mem_points                    = sizeof(Float64) * 3 * size(points, 2)
        b_mem_simps                     = sizeof(Int32) * size(simplices, 1) * size(simplices, 2)
        b_mem_cnt                       = sizeof(UInt16)  * n_samples_x * n_samples_y
        b_mem_misc                      = 16 * 1024^2 * n_threads # 16 MiB per thread
        b_mem_misc += n_samples_x * sizeof(Float64)
        b_mem_misc += n_simplices * sizeof(UInt16)
        b_mem_per_out_var_and_timestep  = sizeof(Float64) * n_samples_x * n_samples_y
        b_mem_per_in_var_and_timestep   = sizeof(Float64) * n_simplices

        n_in_vars                       = length(in_var_names)

        # u, v for 3d; η (eta) for surface; d for seafloor
        out_vars_dyn = 
            if     z_range == z_all;     [(id=1, name="u"), (id=2, name="v")]
            elseif z_range == z_surface; [(id=1, name="eta")]
            else                         [(id=1, name="d")] end

        # b for seafloor
        out_vars_stat = if z_range == z_floor; [(id=1, name="b")] else [] end

        n_out_vars_dyn              = length(out_vars_dyn)
        n_out_vars_stat             = length(out_vars_stat)

        n_timesteps                 = t_end - t_begin + 1
        n_timesteps_per_iteration   = ((mem_limit - b_mem_points - b_mem_simps - b_mem_cnt - b_mem_misc - n_out_vars_stat * b_mem_per_out_var_and_timestep) 
            ÷ (n_in_vars * b_mem_per_in_var_and_timestep + n_out_vars_dyn * b_mem_per_out_var_and_timestep))
        n_timesteps_per_iteration   = max(1, n_timesteps_per_iteration) # Process at least 1 timestep

        n_iterations                = ceil(Int, n_timesteps / n_timesteps_per_iteration)
        
        n_total_problems            = n_timesteps * n_simplices

        s_simplex_name              = n_dims == 3 ? "tetrahedron" : "triangle"
        s_simplex_name_plural       = n_dims == 3 ? "tetrahedra" : "triangles"

        #============================================#
        # Put simplices into bins to avoid thread
        # synchronization & lock contention later on.
        #============================================#

        n_threads = nthreads()
        println("Binning $n_simplices $s_simplex_name_plural into $(n_threads)+$(n_threads - 1)+1 buckets...")
        println("Load balancing using the $(load_balancer == naive ? "naive" : load_balancer == count ? "count" : "workload") balancer.")

        l_x_simps = zeros(Float64, n_samples_x)
        m43_simp_points = Array{Float64, 2}(undef, (n_simplex_points, 3))

        if load_balancer != naive
            #============================================#
            # Ignore this branch, only the naive balancer
            # performs as desired.
            #============================================#

            for tet_id ∈ 1:n_simplices
                simp = (@view simplices[:, tet_id])
                for i ∈ 1:n_simplex_points, dim ∈ 1:3
                    m43_simp_points[i, dim] = points[dim, simp[i]]
                end

                if z_range != z_all
                    if minimum(m43_simp_points[:, 3]) > i_domain_z[2] || maximum(m43_simp_points[:, 3]) < i_domain_z[1]
                        continue
                    end
                end

                x_min = minimum(m43_simp_points[:, 1])
                x_max = maximum(m43_simp_points[:, 1])
                y_min = minimum(m43_simp_points[:, 2])
                y_max = maximum(m43_simp_points[:, 2])
                z_min = minimum(m43_simp_points[:, 3])
                z_max = maximum(m43_simp_points[:, 3])

                idx_x_min = floor(Int, (x_min - i_domain_x[1]) / sampling_rate[1]) + 1
                idx_x_max = floor(Int, (x_max - i_domain_x[1]) / sampling_rate[1]) + 1
                idx_y_min = floor(Int, (y_min - i_domain_y[1]) / sampling_rate[2]) + 1
                idx_y_max = floor(Int, (y_max - i_domain_y[1]) / sampling_rate[2]) + 1
                idx_z_min = floor(Int, (z_min - i_domain_z[1]) / sampling_rate[3]) + 1
                idx_z_max = floor(Int, (z_max - i_domain_z[1]) / sampling_rate[3]) + 1

                idx_x_min = max(1, min(idx_x_min, n_samples_x))
                idx_x_max = max(1, min(idx_x_max, n_samples_x))

                if load_balancer == workload
                    aabb_volume = (idx_x_min - idx_x_min + 1) * (idx_y_max - y_min + 1) * (idx_z_max - idx_z_min + 1)
                    # Obtained through linear regression on benchmark data
                    estimated_workload_factor = lb_params[1] + lb_params[2] * aabb_volume
                    l_x_simps[idx_x_min:idx_x_max] .+= estimated_workload_factor / 1e6 # prevent numerical inaccuracies by downscaling the values
                elseif load_balancer == count
                    l_x_simps[idx_x_min:idx_x_max] .+= 1e-3
                end
            end # for
        else
            # The naive balancer creates buckets of equal size by setting each column to a constant load factor
            l_x_simps .= 1.
        end # if load_balancer != naive

        # Now that each column has a load factor, size each bin such that the load is distributed evenly.
        l_bin_start_idxs = Array{Int, 1}(undef, n_threads)
        n_opt_simps_per_p_bin = sum(l_x_simps) / n_threads
        l_bin_start_idxs[1] = 1

        # Start from the left bin and add columns to it until its rounded load factor is closest to the optimal value (1/n_threads*total_load)
        # Then continue with next bin
        next_x_idx = 1
        for bin_id ∈ 1:(n_threads-1)
            current_simp_count = 0.
            new_simp_count = 0.
            while (n_opt_simps_per_p_bin - current_simp_count > new_simp_count - n_opt_simps_per_p_bin) && next_x_idx <= n_samples_x
                current_simp_count = new_simp_count
                new_simp_count = current_simp_count + l_x_simps[next_x_idx]
                next_x_idx += 1
            end

            l_bin_start_idxs[bin_id + 1] = min(next_x_idx, n_samples_x)
        end

        #============================================#
        # Perform binning
        #============================================#

        # For each tetrahedron, assign a bin with an ID ∈ 1:2*n_threads based on its location
        l_bin_ids = Array{UInt16, 1}(undef, n_simplices)

        Threads.@threads for thread_id ∈ 1:nthreads()
            thread_range_min = floor(INDEX_TYPE, (thread_id - 1) * n_simplices_per_thread) + 1
            thread_range_max = thread_id == n_threads ? n_simplices : floor(INDEX_TYPE, thread_id * n_simplices_per_thread)

            m43_tet_points = Array{Float64, 2}(undef, (n_simplex_points, 3))
            m32_tet_aabb   = Array{Float64, 2}(undef, (3, 2))

            for tet_id ∈ thread_range_min:thread_range_max
                tetrahedron = (@view simplices[:, tet_id])
                for i ∈ 1:n_simplex_points, dim ∈ 1:3
                    m43_tet_points[i, dim] = points[dim, tetrahedron[i]]
                end

                # Simplex is not in the z-domain (only possible if z_range ∈ [z_floor, z_surface]), ignore it.
                if z_range != z_all
                    coord_z_min = minimum(m43_tet_points[:, 3])
                    coord_z_max = maximum(m43_tet_points[:, 3])

                    if coord_z_min > i_domain_z[2] || coord_z_max < i_domain_z[1]
                        l_bin_ids[tet_id] = n_threads * 2 + 1 # Unused bin, will be ignored later
                        continue
                    end
                else # This would be were tetrahedra in the volume mesh that belong to the ground are filtered out. For now: cutoff below -2000m
                    coord_z_min = minimum(m43_tet_points[:, 3])
                    coord_z_max = maximum(m43_tet_points[:, 3])
                    if coord_z_max <= 0.00001
                        l_bin_ids[tet_id] = n_threads * 2 + 1 # Unused bin, will be ignored later
                        continue
                    end
                end

                coord_x_min = minimum(m43_tet_points[:, 1])
                coord_x_max = maximum(m43_tet_points[:, 1])

                idx_x_min = floor(Int,(coord_x_min - i_domain_x[1]) / sampling_rate[1]) + 1
                idx_x_max =  ceil(Int,(coord_x_max - i_domain_x[1]) / sampling_rate[1])

                bin_id_l = findlast(bin_start_idx -> bin_start_idx ≤ idx_x_min, l_bin_start_idxs)
                bin_id_r = findlast(bin_start_idx -> idx_x_max ≥ bin_start_idx, l_bin_start_idxs)
                if bin_id_r > n_threads; bin_id_r = n_threads; end # The right border of the last bin is inclusive while the others are not.

                # n_bins := 2*n_threads
                # bin i ∈ 1:n_bins contains the tets completely inside the column of bin i
                # bin j ∈ n_bins+1:2*n_bins-1 contains tets that are overlapping across bin i and i+1
                # bin k = 2*n_bins contains tets overlapping over more than 2 adjascent columns
                if bin_id_l == bin_id_r
                    l_bin_ids[tet_id] = bin_id_l
                elseif bin_id_l == bin_id_r - 1
                    l_bin_ids[tet_id] = n_threads + bin_id_l
                else
                    l_bin_ids[tet_id] = 2*n_threads
                end
            end # for tet_id
        end # for thread_id

        #============================================#
        # Calculate number of simps in each bin,
        # print statistics
        #============================================#

        # The number of tetrahedra in each respective bin (only the first 2*n_threads bins have been filled)
        l_bin_counts = zeros(INDEX_TYPE, 2*n_threads + 1)

        for bin_id ∈ l_bin_ids
            l_bin_counts[bin_id] += 1
        end

        for bin_id ∈ 1:length(l_bin_counts) - 1
            println(l_bin_counts[bin_id], " ", s_simplex_name_plural, " in bin ", bin_id)
        end

        if z_range != z_all
            println(l_bin_counts[length(l_bin_counts)], " ", s_simplex_name_plural, " are ignored due to z_range filter.")
        end

        println("Done.")

        #============================================#
        # Set up NetCDF output file
        #============================================#

        if !lb_autotune # LB-autotune does not write outputs
            if create_file 
                Main.NC.get_or_create_netcdf(out_filename; create=true,
                                                y_vals=[i_domain_y[1] + i * sampling_rate[2] for i ∈ 1:n_samples_y],
                                                x_vals=[i_domain_x[1] + i * sampling_rate[1] for i ∈ 1:n_samples_x],
                                                t_vals=times[t_begin:t_end])
            else
                Main.NC.get_or_create_netcdf(out_filename)
            end
        end

        #============================================#
        # Setup output grids in memory
        #============================================#

        # These are the grids that will be written to the NetCDF output later on.
        l_dyn_output_grids = Array{Array{Float64, 3}, 1}(undef, n_out_vars_dyn)
        for var ∈ out_vars_dyn
            l_dyn_output_grids[var.id] = Array{Float64, 3}(undef, (n_samples_x, n_samples_y, n_timesteps_per_iteration))
        end

        l_stat_output_grids = Array{Array{Float64, 2}, 1}(undef, n_out_vars_stat)
        for var ∈ out_vars_stat
            l_stat_output_grids[var.id] = Array{Float64, 2}(undef, (n_samples_x, n_samples_y))
        end

        # This matrix counts the total number of samples in each output grid cell
        # (read: the total number of sampled points in z-direction for each xy-index of the output grids)
        # This stays constant across all vars and timesteps because the grid geometry is constant.
        myx_output_sample_counts = Array{UInt16, 2}(undef, (n_samples_x, n_samples_y))

        #============================================#
        # Sample tetrahedra and write to output
        #============================================#

        println("Processing timesteps ", t_begin-1, "-", t_end-1, " of 0-", length(times)-1, ".")
        println("RAM limit: ", Main.Util.human_readable_size(mem_limit), "; Can work on ", 
                n_timesteps_per_iteration, " timesteps at once (therefore needing ", 
                n_iterations , " iterations).")

        # In each iteration, process ALL tetrahedra for ALL variables for as many timesteps as possible*
        # *possible means that the maximum memory footprint set above cannot be exceeded. This limits the number of 
        # timesteps that can be stored in l_iter_output_grids. 
        for iteration ∈ 1:n_iterations
            t_start = t_begin + (iteration - 1) * n_timesteps_per_iteration
            n_times = min(n_timesteps - (iteration - 1) * n_timesteps_per_iteration, n_timesteps_per_iteration)

            #============================================#
            # Prefetch the values of all vars in all
            # timesteps for the current iteration.
            # This eliminates random accesses later on
            # and D R A S T I C A L L Y improves runtime.
            #============================================#

            prefetched_vars = Array{Float64, 3}(undef, (n_simplices, n_times, n_in_vars))
            @profile begin
                for var_id ∈ 1:n_in_vars, time ∈ t_start:t_start + n_times - 1
                    dest_offset = (var_id-1)*n_times*n_simplices + (time-t_start)*n_simplices + 1
                    copyto!(prefetched_vars, dest_offset, in_vars[time, var_id], 1, n_simplices)
                end
            end # profile

            # Reset the output values to 0. Only do so for the timesteps actually processed in this iteration.
            for grid ∈ l_dyn_output_grids
                grid[:,:, 1:n_times] .= 0
            end

            # Static vars only need to be read once. Do that in the first iteration.
            if iteration == 1
                for grid ∈ l_stat_output_grids
                    grid .= 0
                end
            end

            myx_output_sample_counts .= 0

            #============================================#
            # Binning
            #============================================#
            # First, process each of the n_threads normal
            # bins. (P-partitioning)
            #
            # Then, process the n_threads-1 bins on the
            # borders between normal bins. (B-partitioning)
            #
            # Lastly, process the 1 bin of tetrahedra
            # that overlapped more than 2 adjacent bins.
            # (R-partitioning)
            #============================================#

            # Allocate structures needed for the load balancing autotune routine (if enabled)
            l_tet_autotune_values :: Union{Nothing, Array{Int, 2}} = nothing
            if lb_autotune
                l_tet_autotune_values = fill(-1, (n_simplices, 2))
            end

            # P-partitioning
            Threads.@threads for thread_id ∈ 1:n_threads
                iterate_tets!(
                    thread_id,         l_bin_ids,    l_bin_counts,
                    n_simplices,       simplices,    points, 
                    prefetched_vars,   out_vars_dyn, out_vars_stat,
                    i_domain_x,        i_domain_y,   i_domain_z,
                    n_samples_x,       n_samples_y,  n_samples_z,
                    sampling_rate,     t_start,      n_times,
                    l_dyn_output_grids, l_stat_output_grids, myx_output_sample_counts, n_simplex_points, z_range,
                    print_progress = (thread_id == n_threads ÷ 2), lb_autotune=l_tet_autotune_values)
            end # for thread_id

            println("Processing $s_simplex_name_plural on bucket borders...")

            # B-partitioning
            Threads.@threads for thread_id ∈ 1:n_threads-1
                iterate_tets!(
                    n_threads + thread_id, l_bin_ids,    l_bin_counts,
                    n_simplices,           simplices,    points, 
                    prefetched_vars,       out_vars_dyn, out_vars_stat,
                    i_domain_x,            i_domain_y,   i_domain_z,
                    n_samples_x,           n_samples_y,  n_samples_z,
                    sampling_rate,         t_start,      n_times,
                    l_dyn_output_grids, l_stat_output_grids, myx_output_sample_counts, n_simplex_points, z_range, 
                    print_progress = (thread_id == n_threads ÷ 2), lb_autotune=l_tet_autotune_values)
            end # for thread_id

            # R-partitioning
            iterate_tets!(
                2 * n_threads,     l_bin_ids,    l_bin_counts,
                n_simplices,       simplices,    points, 
                prefetched_vars,   out_vars_dyn, out_vars_stat,
                i_domain_x,        i_domain_y,   i_domain_z,
                n_samples_x,       n_samples_y,  n_samples_z,
                sampling_rate,     t_start,      n_times,
                l_dyn_output_grids, l_stat_output_grids, myx_output_sample_counts, n_simplex_points, z_range, 
                print_progress = true, lb_autotune=l_tet_autotune_values)

            println("Done.")

            # Collect the AABB volumes and runtimes for each simplex and calculate the workload-LB params through linear regression.
            # Currently does not deliver desirable results. Do not use.
            if lb_autotune
                open("linload.csv", "w+") do fp
                    println("aabbVol;aabbTime")
                    for i in 1:size(l_tet_autotune_values, 1)
                        println(fp, l_tet_autotune_values[i, 1], ';', l_tet_autotune_values[i, 2])
                    end
                end
                lb_params = [ones(n_simplices) l_tet_autotune_values[:, 1]] \ l_tet_autotune_values[:, 2]
                println("==========================================")
                println("  workload-LB params: a = $(lb_params[1]), b = $(lb_params[2])")
                println("==========================================")

                write("sampler_lb_params.txt", string(lb_params[1]), ';', string(lb_params[2]))

                println()
                println("LB-params written to file, will use them automatically from now on when --load-balancer=workload is specified.")
                exit(0) # LB-autotune complete, kill program
            end

            #============================================#
            # Write NetCDF outputs for this iteration
            #============================================#

            println("Writing outputs for timesteps $(t_start-1)-$(t_start+n_times-2)...")
            
            start = [1, 1, (t_start-t_begin)+1]
            count = [-1, -1, n_times]

            for var ∈ out_vars_dyn
                ncwrite(l_dyn_output_grids[var.id], out_filename, var.name, start=start, count=count)
            end

            if iteration == 1 # Only write static outputs once
                for var ∈ out_vars_stat
                    if var.name == "b" # Water height correction
                        l_stat_output_grids[var.id] .-= water_height
                    end
                    ncwrite(l_stat_output_grids[var.id], out_filename, var.name, start=[1, 1], count=[-1, -1])
                end
            end

            println("Done.")
        end # for iteration

    end # function rasterize

    function iterate_tets!(
        bin_id,           l_bin_ids,    l_bin_counts, 
        n_tetrahedra,     l_simplices, l_points, 
        in_vars,          out_vars_dyn, out_vars_stat,
        i_domain_x,       i_domain_y,   i_domain_z,
        n_samples_x,      n_samples_y,  n_samples_z,
        v_sampling_rate,  t_start,      n_times,
        l_dyn_grids,      l_stat_grids, myx_grid_sample_counts, 
        n_simplex_points, z_range;      print_progress=false, lb_autotune :: Union{Nothing, Array{Int, 2}} = nothing)

        #============================================#
        # Pre-allocate memory used throughout
        # rasterization by this thread
        #============================================#

        m43_simp_points      = Array{Float64, 2}(undef, (n_simplex_points, 3))
        m32_tet_aabb        = Array{Float64, 2}(undef, (3, 2))
        l_face_points       = Array{SubArray, 1}(undef, n_simplex_points)
        l_tet_faces         = Array{HesseNormalForm, 1}(undef, n_simplex_points)
        
        v_p                 = Array{Float64, 1}(undef, 3)
        v_n                 = Array{Float64, 1}(undef, 3)
        simp_sample_counts   = Array{UInt16, 2}(undef, (64, 64))
        
        d_start             = now()
        d_last_printed      = d_start
        print_interval      = Second(2)
        n_bin_rasterized    = 0
        n_bin_total         = l_bin_counts[bin_id]

        for tet_id ∈ 1:n_tetrahedra
            # Only process tetrahedra in thread's own bin
            if l_bin_ids[tet_id] != bin_id
                continue
            end

            if print_progress
                d_now = now()

                if (d_now - d_last_printed) >= print_interval && n_bin_rasterized > 0
                    d_last_printed = d_now

                    etr  = (d_now-d_start) * (n_bin_total - n_bin_rasterized) ÷ (n_bin_rasterized)
                    hh   = floor(etr, Dates.Hour)
                    etr -= floor(hh, Dates.Millisecond)
                    mm   = floor(etr, Dates.Minute)
                    etr -= floor(mm, Dates.Millisecond)
                    ss   = floor(etr, Dates.Second)

                    @printf("Working on simplex %d of %d in timesteps %d-%d (%2.2f%% done). ETR: %02d:%02d:%02d\n",
                            n_bin_rasterized, n_bin_total,
                            t_start - 1,           # Tools like ParaView work 0-indexed. Avoid confusion by outputting 0-indexed here.
                            t_start + n_times - 2, # Tools like ParaView work 0-indexed. Avoid confusion by outputting 0-indexed here.
                            n_bin_rasterized/n_bin_total*100,
                            hh.value, mm.value, ss.value)
                end
            end # if print_progress

            timeit_t_start = time_ns() # LB-autotune

            # Read simplex points. m43 is only filled with 3x3 values for triangles
            l_simplex = (@view l_simplices[:, tet_id])
            for i ∈ 1:n_simplex_points, dim ∈ 1:3
                m43_simp_points[i, dim] = l_points[dim, l_simplex[i]]
            end

            for dim ∈ 1:3
                m32_tet_aabb[dim, 1] = minimum(m43_simp_points[:, dim])
                m32_tet_aabb[dim, 2] = maximum(m43_simp_points[:, dim])
            end

            #============================================#
            # A note on global & local indices
            #============================================#
            # In the following, idx_g_ will refer to an
            # index of the global sampling domain while
            # idx_l_ refers to an index in the simp-local
            # domain.
            #============================================#

            idx_g_x_min = floor(Int32, (m32_tet_aabb[1,1] - i_domain_x[1]) / v_sampling_rate[1]) + 1
            idx_g_y_min = floor(Int32, (m32_tet_aabb[2,1] - i_domain_y[1]) / v_sampling_rate[2]) + 1
            idx_g_z_min = floor(Int32, (m32_tet_aabb[3,1] - i_domain_z[1]) / v_sampling_rate[3]) + 1
     
            idx_g_x_max =  ceil(Int32, (m32_tet_aabb[1,2] - i_domain_x[1]) / v_sampling_rate[1])
            idx_g_y_max =  ceil(Int32, (m32_tet_aabb[2,2] - i_domain_y[1]) / v_sampling_rate[2])
            idx_g_z_max = floor(Int32, (m32_tet_aabb[3,2] - i_domain_z[1]) / v_sampling_rate[3]) + 1

            n_current_cells_x = idx_g_x_max - idx_g_x_min + 1
            n_current_cells_y = idx_g_y_max - idx_g_y_min + 1
            n_current_cells_z = idx_g_z_max - idx_g_z_min + 1

            # If the current simp has a bounding box larger than anticipated, enlarge the sample counts array accordingly
            # By allocating double the size of the bounding box, future simps that are just a bit larger than the current one
            # will not trigger a resize operation.
            if size(simp_sample_counts, 1) < n_current_cells_y || size(simp_sample_counts, 2) < n_current_cells_x
                simp_sample_counts = Array{UInt16, 2}(undef, (n_current_cells_y * 2, n_current_cells_x * 2))
            end

            # For the bounding box of the current tet, reset all sample counts to 0
            simp_sample_counts[1:n_current_cells_y, 1:n_current_cells_x] .= 0

            for i ∈ 1:n_simplex_points
                l_face_points[i] = @view m43_simp_points[i,:]
            end

            for i ∈ 1:n_simplex_points
                # Exclude point i from the tetrahedron and consider the plane defined by the other 3 points
                # If and only if (x,y,z) lies on the "positive" side (in normal direction) of all 4 planes, 
                # it is inside the tetrahedron.
                # Analogous for triangles but the fourth point is placed straight above the third one which
                # allows the usage of the same point-in-simplex code afterwards.
                v_p_excl = l_face_points[i]
                v_p1 =     l_face_points[( i    % n_simplex_points) + 1]
                v_p2 =     l_face_points[((i+1) % n_simplex_points) + 1]

                v_p3 = 
                    if n_simplex_points == 4; l_face_points[((i+2) % n_simplex_points) + 1]
                    elseif n_simplex_points == 3; [v_p1[1], v_p1[2], v_p1[3] + 1.]
                    else throw(DimensionMismatch("Only simplices with 3 or 4 points (read: triangles / tetrahedra) are supported."))
                    end
            
                # Calculate Hesse normal form of the plane defined by p1, p2, p3
                # Normal vector
                cross3!(v_p2-v_p1, v_p3-v_p1, v_n)
                n_simplex_points == 3 && @assert (v_n[3] == 0.) "Normal vector of 2D line has non-zero z component!"
                normalize!(v_n)
            
                # n1x1 + n2x2 + n3x3 + d = 0; solve for d
                d = -dot3(v_p1, v_n)
                dist_p_excl = dot3(v_p_excl, v_n) + d

                # Ensure that p_excl is considered to be on the "positive" side of the plane
                if dist_p_excl < 0
                    v_n = -v_n
                    d = -d
                end

                l_tet_faces[i] = HesseNormalForm((v_n[1], v_n[2], v_n[3]), d)
            end # for i

            @inline function is_in_simplex()
                accepted = true
                # For each tetrahedron face, filter out cells that are not in the tetrahedron
                for i ∈ 1:n_simplex_points
                    face = l_tet_faces[i]
                    # p_excl is ALWAYS on the "positive" side of the plane. 
                    # So, if p_excl and p lie on the same side, p is inside the tetrahedron
                    dist_p = dot3(v_p, face.n0) + face.d
                
                    # The point lies approx. ON the face, reject in the tet on the lexicographically
                    # more negative side of the face. This ensures that for two neighboring tetrahedra, the point is
                    # only accepted for one of them.
                    # Def.: lexicographical order: 
                    #       a ≤ b ⟺ (a.x ≤ b.x) ∨ (a.x = b.x ∧ a.y ≤ b.y) ∨ (a.x = b.x ∧ a.y = b.y ∧ a.z ≤ b.z)
                    # This is a total ordering on all vectors in R^3 and will therefore always be larger for the "positive"
                    # normal vector in contrast to its negated counterpart.
                    if (abs(dist_p) <= EDGE_TOL) 
                        is_positive_side = 
                            (face.n0[1] > 0.) || 
                            (face.n0[1] == 0. && face.n0[2] > 0.) || 
                            (face.n0[1] == 0. && face.n0[2] == 0. && face.n0[3] > 0.)
                    
                        # Reject if in tetrahedron on face's lex. negative side
                        if !is_positive_side
                            accepted = false
                            break
                        end
                    end
                
                    # If signs of distance vectors do not match, p and p_excl are on different sides of the plane (reject)
                    # As stated above, the normal vector n is chosen so that p_excl is on the POSITIVE side of the plane
                    # Therefore, we only need to check if dist_p is negative to reject p.
                    if dist_p < EDGE_TOL
                        accepted = false
                        break
                    end
                end # for i

                return accepted
            end

            #============================================#
            # Sample each cuboid that is overlapping the 
            # current simplex's AABB.
            # If its CENTER POINT lies within the simp,
            # increase the simp's sample count for the
            # cuboid's (x, y) index by 1.
            #============================================#

            for idx_g_x ∈ idx_g_x_min:idx_g_x_max
                v_p[1] = (idx_g_x-1)*v_sampling_rate[1] + i_domain_x[1]
                for idx_g_y ∈ idx_g_y_min:idx_g_y_max
                    v_p[2] = (idx_g_y-1)*v_sampling_rate[2] + i_domain_y[1]

                    if n_simplex_points == 4 # processing tets, increment counter for intersected cells
                        for idx_g_z ∈ idx_g_z_min:idx_g_z_max
                            v_p[3] = (idx_g_z-1)*v_sampling_rate[3] + i_domain_z[1]
                            if is_in_simplex() # accept cell and increment its 2D cell counter for the current tetrahedron by 1
                                idx_l_x = idx_g_x - idx_g_x_min + 1
                                idx_l_y = idx_g_y - idx_g_y_min + 1
                                simp_sample_counts[idx_l_y, idx_l_x] += 1
                            end
                        end # for z
                    else # processing triangles, immediately rasterize vars
                        if is_in_simplex()
                            idx_l_x = idx_g_x - idx_g_x_min + 1
                            idx_l_y = idx_g_y - idx_g_y_min + 1

                            if z_range == z_floor
                                # Only calculate bathymetry height in the first timestep since it is static
                                if t_start == 1
                                    v_λ = bary_coords_2d(v_p, m43_simp_points)

                                    # Calculate bathymetry height from triangle's z-coords and interpolate between them
                                    b = sum(m43_simp_points[1:3,3] .* v_λ) / 3.
                                    l_stat_grids[1][idx_g_x, idx_g_y] = b
                                end

                                for t ∈ 1:n_times
                                    l_dyn_grids[1][idx_g_x, idx_g_y, t] = in_vars[tet_id, t, 1]
                                end
                            elseif z_range == z_surface
                                for t ∈ 1:n_times
                                    η = in_vars[tet_id, t, 1]
                                    l_dyn_grids[1][idx_g_x, idx_g_y, t] = η
                                end
                            else # Processing triangles but neither on floor nor on surface. No use for that at the moment
                                throw(ErrorException("Not implemented."))
                            end

                            simp_sample_counts[idx_l_y, idx_l_x] = 1 # 2D surface is only sampled once in z-direction
                        end
                    end
                end # for y
            end # for x

            #============================================#
            # All cuboids overlapping the tet's domain
            # have been sampled, 
            # now calculate the running average of each
            # (x, y) output grid cell with the newly gained
            # samples.
            # The weight of each tetrahedron in the
            # average is the number of samples of the grid
            # cell lying within the tet divided by the total
            # number of samples in the cell.
            #============================================#

            if z_range == z_all # Build running average of variables in volume grid
                for idx_l_x ∈ 1:n_current_cells_x, idx_l_y ∈ 1:n_current_cells_y
                    num_samples = simp_sample_counts[idx_l_y, idx_l_x]
                    if num_samples > 0
                        idx_g_x = idx_l_x + idx_g_x_min - 1; idx_g_y = idx_l_y + idx_g_y_min - 1
                        total_samples = (myx_grid_sample_counts[idx_g_x, idx_g_y] += num_samples) # x, y flip intentional.
                        for var ∈ out_vars_dyn
                            for t ∈ 1:n_times
                                avg = l_dyn_grids[var.id][idx_g_x, idx_g_y, t] # x, y flip intentional
                                l_dyn_grids[var.id][idx_g_x, idx_g_y, t] = avg + (in_vars[tet_id, t, var.id]-avg)*(num_samples/total_samples)
                            end
                        end
                    end
                end
            end # if z_range

            n_bin_rasterized += 1
            if !isnothing(lb_autotune)
                timeit_time = time_ns() - timeit_t_start

                lb_autotune[tet_id, 1] = n_current_cells_x * n_current_cells_y * n_current_cells_z
                lb_autotune[tet_id, 2] = timeit_time

                if n_simplex_points == 4; lb_autotune[tet_id, 1] *= n_current_cells_z; end
            end
        end # for tet_id
    end

    """
    In-place cross product to avoid heaps of heap allocations
    """
    @inline function cross3!(a, b, ret)
        ret[1] = (a[2]*b[3] - a[3]*b[2])
        ret[2] = (a[3]*b[1] - a[1]*b[3])
        ret[3] = (a[1]*b[2] - a[2]*b[1])
    end

    """
    Explicit dot product for 3D vectors to avoid huge call stacks
    """
    @inline function dot3(a, b)
        return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
    end

    """
    Calculate barycentric coordinates of `p` in the triangle defined by `ps`.
    """
    @inline function bary_coords_2d(p :: Array{Float64, 1}, ps :: Array{Float64, 2}) :: NTuple{3, Float64}
        @assert length(p) >= 2
        @assert size(ps, 1) >= 3
        @assert size(ps, 2) >= 2

        T11 = ps[1, 1] - ps[3, 1]
        T12 = ps[2, 1] - ps[3, 1]
        T21 = ps[1, 2] - ps[3, 2]
        T22 = ps[2, 2] - ps[3, 2]

        det = T11 * T22 - T12 * T21

        λ1 = (T22*(p[1]-ps[3,1]) - T12*(p[2]-ps[3,2])) / det
        λ2 = (T11*(p[2]-ps[3,2]) - T21*(p[1]-ps[3,1])) / det
        λ3 = 1. - λ1 - λ2

        return (λ1, λ2, λ3)
    end
end # module Rasterization
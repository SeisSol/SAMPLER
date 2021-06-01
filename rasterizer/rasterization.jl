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
    using Main.XDMF
    using Main.Args
    using Infiltrator

    export rasterize, ZRange

    @enum ZRange z_all=-1 z_floor=0 z_surface=1

    const EDGE_TOL = .0      # Tolerance for whether to consider points that are close to a simplex edge/face to be on that edge/face
    const Z_RANGE_TOL = .001 # Tolerance for points close to the z_range-limit. Points over the limit by less than Z_RANGE_TOL are not discarded

    # Aliases for common array indices
    const X = 1
    const Y = 2
    const Z = 3
    const MIN = 1
    const MAX = 2

    struct HesseNormalForm
        n0  :: NTuple{3, Float64}
        d   :: Float64
    end

    struct RasterizationContext
        xdmf                        :: XDMF.XDMFFile
        simplices                   :: AbstractArray{Integer, 2}
        points                      :: AbstractArray{AbstractFloat, 2}

        z_range                     :: ZRange
        domain                      :: NTuple{3, NTuple{2, Float64}} # (x, y, z) with (min, max) respectively
        samples                     :: NTuple{3, UInt}
        sampling_rate               :: NTuple{3, Float64}
        dyn_var_mapping             :: Main.Args.VarMapping
        stat_var_mapping            :: Main.Args.VarMapping
        in_vars_dyn                 :: Array{String, 1}
        out_vars_stat               :: Array{String, 1}
        out_vars_dyn                :: Array{String, 1}

        n_threads                   :: Int
        n_simplex_points            :: Int
        n_dims                      :: Int
        n_simplices                 :: Int
        n_in_vars                   :: Int
        n_out_vars_dyn              :: Int
        n_out_vars_stat             :: Int
        n_timesteps                 :: Int
        n_timesteps_per_iteration   :: Int
        n_iterations                :: Int

        simplices_per_thread        :: Float64

        tanioka                     :: Bool
        water_height                :: Float64

        t_begin                     :: UInt
        t_end                       :: UInt
        times                       :: AbstractArray{Float64, 1}

        out_filename                :: AbstractString
    end

    struct IterationBuffers
        bin_ids
        bin_counts
        out_grids_dyn               :: Dict{String, Array{Float64, 3}}
        out_grids_stat              :: Dict{String, Array{Float64, 2}}
        out_sample_counts           :: Array{UInt16, 2}
        prefetched_vars             :: Dict{String, Array{Float64, 2}}
    end

    @inline function simplex_name(ctx, plural=true)
        return ctx.n_dims == 3 ? (plural ? "simplices" : "simplex") : (plural ? "triangles" : "triangle") 
    end

    function RasterizationContext(xdmf, sampling_rate, mem_limit, dyn_var_mapping, stat_var_mapping, z_range, tanioka, water_height, t_begin, t_end, times, out_filename)
        simplices, points               = grid_of(xdmf)
        domain                          = Tuple( (minimum(points[i, :]), maximum(points[i, :]))                     for i in X:Z)
        samples                         = Tuple( ceil(UInt, (domain[i][MAX] - domain[i][MIN]) / sampling_rate[i])   for i in X:Z)

        #============================================#
        # Calculate problem size & hardware limitations
        #============================================#

        n_threads                       = nthreads()

        n_simplex_points                = size(simplices, 1)
        n_dims                          = n_simplex_points - 1
        n_simplices                     = size(simplices, 2)
        simplices_per_thread            = n_simplices / n_threads # float on purpose

        b_mem_points                    = sizeof(Float64) * 3 * size(points, 2)
        b_mem_simps                     = sizeof(Int32) * size(simplices, 1) * size(simplices, 2)
        b_mem_cnt                       = sizeof(UInt16)  * samples[X] * samples[Y]
        b_mem_misc                      = 16 * 1024^2 * n_threads # 16 MiB per thread
        b_mem_misc                     += samples[X] * sizeof(Float64)
        b_mem_misc                     += n_simplices * sizeof(UInt16)
        b_mem_per_out_var_and_timestep  = sizeof(Float64) * samples[X] * samples[Y]
        b_mem_per_in_var_and_timestep   = sizeof(Float64) * n_simplices

        n_in_vars                       = length(dyn_var_mapping)
        # static vars are currently only bathymetry and thus extracted from grid geometry. Thus, no in_vars_stat.
        in_vars_dyn                     = collect(keys(dyn_var_mapping))

        out_vars_stat                   = z_range != z_all ? collect(values(stat_var_mapping)) : []
        out_vars_dyn                    = collect(values(dyn_var_mapping))

        n_out_vars_dyn                  = length(dyn_var_mapping)
        n_out_vars_stat                 = length(out_vars_stat)

        n_timesteps                     = t_end - t_begin + 1
        n_timesteps_per_iteration       = ((mem_limit - b_mem_points - b_mem_simps - b_mem_cnt - b_mem_misc - n_out_vars_stat * b_mem_per_out_var_and_timestep)
                                        ÷  (n_in_vars * b_mem_per_in_var_and_timestep + n_out_vars_dyn * b_mem_per_out_var_and_timestep))
        n_timesteps_per_iteration       = max(1, n_timesteps_per_iteration) # Process at least 1 timestep

        n_iterations                    = ceil(Int, n_timesteps / n_timesteps_per_iteration)

        return RasterizationContext(xdmf, simplices, points, 
                                    z_range, domain, samples, sampling_rate,
                                    dyn_var_mapping, stat_var_mapping, 
                                    in_vars_dyn, out_vars_stat, out_vars_dyn,
                                    n_threads, n_simplex_points, 
                                    n_dims, n_simplices, 
                                    n_in_vars, n_out_vars_dyn, 
                                    n_out_vars_stat, n_timesteps,
                                    n_timesteps_per_iteration, n_iterations, 
                                    simplices_per_thread,
                                    tanioka, water_height, t_begin, t_end, times,
                                    out_filename)
    end

    function IterationBuffers(bin_ids, bin_counts, ctx :: RasterizationContext)
        # These are the grids that will be written to the NetCDF output later on.
        out_grids_dyn = Dict{String, Array{Float64, 3}}()
        for var_name ∈ ctx.out_vars_dyn
            out_grids_dyn[var_name] = Array{Float64, 3}(undef, (ctx.samples[X:Y]..., ctx.n_timesteps_per_iteration))
        end

        out_grids_stat = Dict{String, Array{Float64, 2}}()
        for var_name ∈ ctx.out_vars_stat
            out_grids_stat[var_name] = Array{Float64, 2}(undef, ctx.samples[X:Y])
        end

        # This matrix counts the total number of samples in each output grid cell
        # (read: the total number of sampled points in z-direction for each xy-index of the output grids)
        # This stays constant across all vars and timesteps because the grid geometry is constant.
        out_sample_counts = Array{UInt16, 2}(undef, ctx.samples[X:Y])

        # Variable prefetch buffers, get rewritten every iteration
        prefetched_vars = Dict{String, Array{Float64, 2}}()
        for var_name ∈ ctx.in_vars_dyn
            prefetched_vars[var_name] = Array{Float64, 2}(undef, (ctx.n_simplices, ctx.n_timesteps_per_iteration))
        end

        return IterationBuffers(bin_ids, bin_counts, out_grids_dyn, out_grids_stat, out_sample_counts, prefetched_vars)
    end

    """
        Rasterize variables stored in an unstructured grid into a regular grid. Write the results to a NetCDF file.

    # Arguments
    - `xdmf`:             The XDMF file handle of the file to be rasterized
    - `in_var_names`:     The names, of the variables in `in_vars` in the same order as those vars
    - `in_out_mapping`:   A dictionary providing the output variable name for each input variable name
    - `times`:            An array of timestamps associated to the timesteps
    - `sampling_rate`:    A tuple (Δx, Δy, Δz) defining the regular grid cell size.
    - `out_filename`:     The path/name of the output NetCDF file
    - `mem_limit`:        The soft memory (in bytes) the software is allowed to utilize (soft limit!)
    - `z_range`:          A filter for discarding simplices above/below a threshold. `z_floor` discards everything above `0-Z_RANGE_TOL`,
                          `z_surface` discards everything below. `z_all` keeps all simplices.
    - `t_begin`:          The first timestep to be rasterized
    - `t_end`:            The last timestep to be rasterized
    - `create_file_vars`: If not empty: create output file containing the var names in the list (already existing file will be overwritten).
                          Else: use existing output file.
    - `water_height`:     The offset to translate the water level to be at `z=0`.
    - `tanioka`:          Whether to apply Tanioka's method to the bathymetry. Requires U, V, W displacements as input.
    """
    function rasterize(
        xdmf                :: XDMFFile,
        dyn_var_mapping     :: Main.Args.VarMapping,
        stat_var_mapping    :: Main.Args.VarMapping,
        times               :: AbstractArray{Float64, 1},
        sampling_rate       :: NTuple{3, Float64},
        out_filename        :: AbstractString,
        mem_limit           :: Int;
        z_range             :: ZRange  = z_all,
        t_begin             :: Integer = 1,
        t_end               :: Integer = length(times),
        create_nc_dyn_vars  :: Array{Main.Args.VarMapping, 1} = Array{Main.Args.VarMapping, 1}(),
        create_nc_stat_vars :: Array{Main.Args.VarMapping, 1} = Array{Main.Args.VarMapping, 1}(),
        water_height        :: Float64 = 0.,
        tanioka             :: Bool    = false)

        ctx = RasterizationContext(xdmf, sampling_rate, mem_limit, dyn_var_mapping, stat_var_mapping, z_range, tanioka, water_height, t_begin, t_end, times, out_filename)

        #============================================#
        # Print rasterization settings
        #============================================#

        if (ctx.z_range == z_all) != (ctx.n_dims == 3)
            error("Z-Range 'all' can only be used with tetrahedra; Z-Ranges 'floor' and 'surface' can only be used with triangles.")
        end

        println("Domain is $(ctx.domain[X]) × $(ctx.domain[Y]) × $(ctx.domain[Z]).")

        if ctx.z_range != z_all
            s_region = ctx.z_range == z_floor ? "floor" : "surface"
            println("Only processing triangles on the $s_region.")
        end

        if ctx.n_dims == 3
            println("Sampling into $(ctx.samples[X]) × $(ctx.samples[Y]) × $(ctx.samples[Z]) cells.")
        else
            println("Sampling into $(ctx.samples[X]) × $(ctx.samples[Y]) cells.")
        end

        #============================================#
        # Set up NetCDF output file
        #============================================#

        if !isempty(create_nc_dyn_vars) || !isempty(create_nc_stat_vars)
            Main.NC.create_netcdf(ctx.out_filename, [ctx.domain[Y][MIN] + i * ctx.sampling_rate[Y] for i ∈ 1:ctx.samples[Y]],
                                  [ctx.domain[X][MIN] + i * ctx.sampling_rate[X] for i ∈ 1:ctx.samples[X]],
                                  ctx.times[ctx.t_begin:ctx.t_end],
                                  create_nc_stat_vars, create_nc_dyn_vars)
        else
            Main.NC.get_netcdf(ctx.out_filename)
        end

        #============================================#
        # Put simplices into bins to avoid thread
        # synchronization & lock contention later on.
        #============================================#

        bin_ids, bin_counts = bin(ctx)

        #============================================#
        # Setup output grids in memory
        #============================================#

        itbuf = IterationBuffers(bin_ids, bin_counts, ctx)

        #============================================#
        # Sample simplices and write to output
        #============================================#

        println("Processing timesteps ", ctx.t_begin-1, "-", ctx.t_end-1, " of 0-", length(times)-1, ".")
        println("RAM limit: ", Main.Util.human_readable_size(mem_limit), "; Can work on ",
                ctx.n_timesteps_per_iteration, " timesteps at once (therefore needing ",
                ctx.n_iterations , " iterations).")

        #============================================#
        # Iterate over and rasterize timesteps
        #============================================#

        iterate(ctx, itbuf)

    end # function rasterize

    function iterate(ctx, itbuf)
        # In each iteration, process ALL simplices for ALL variables for as many timesteps as possible*
        # *possible means that the maximum memory footprint set above cannot be exceeded. This limits the number of
        # timesteps that can be stored in l_iter_output_grids.
        for iteration ∈ 1:ctx.n_iterations
            t_start = ctx.t_begin + (iteration - 1) * ctx.n_timesteps_per_iteration
            n_times = min(ctx.n_timesteps - (iteration - 1) * ctx.n_timesteps_per_iteration, ctx.n_timesteps_per_iteration)
            t_stop = t_start + n_times - 1

            #============================================#
            # Prefetch the values of all vars in all
            # timesteps for the current iteration.
            # This eliminates random accesses later on
            # and D R A S T I C A L L Y improves runtime.
            #============================================#

            var_bufs :: Dict{AbstractString, Array{AbstractArray, 1}} = Main.XDMF.data_of(ctx.xdmf, t_start, t_stop, ctx.in_vars_dyn)

            for var_name ∈ ctx.in_vars_dyn, t ∈ 1:n_times
                copyto!(itbuf.prefetched_vars[var_name][1:ctx.n_simplices, t], var_bufs[var_name][t][1:ctx.n_simplices])
            end

            # Reset the output values to 0. Only do so for the timesteps actually processed in this iteration.
            for grid ∈ values(itbuf.out_grids_dyn)
                grid[:,:, 1:n_times] .= 0
            end

            # Static vars only need to be read once. Do that in the first iteration.
            if iteration == 1
                for grid ∈ values(itbuf.out_grids_stat)
                    grid .= 0
                end
            end

            itbuf.out_sample_counts .= 0

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

            # P-partitioning
            Threads.@threads for thread_id ∈ 1:ctx.n_threads
                iterate_simps!(thread_id, ctx, itbuf, t_start, t_stop, print_progress = (thread_id == ctx.n_threads ÷ 2))
            end # for thread_id

            println("Processing $(simplex_name(ctx)) on bucket borders...")

            # B-partitioning
            Threads.@threads for thread_id ∈ 1:ctx.n_threads-1
                iterate_simps!(ctx.n_threads + thread_id, ctx, itbuf, t_start, t_stop, print_progress = (thread_id == ctx.n_threads ÷ 2))
            end # for thread_id

            # R-partitioning
            iterate_simps!(2 * ctx.n_threads, ctx, itbuf, t_start, t_stop, print_progress = true)

            println("Done.")

            #============================================#
            # Write NetCDF outputs for this iteration
            #============================================#

            println("Writing outputs for timesteps $(t_start-1)-$(t_stop-1)...")

            start = [1, 1, (t_start-ctx.t_begin)+1]
            cnt = [-1, -1, n_times]

            for name ∈ ctx.out_vars_dyn
                ncwrite(itbuf.out_grids_dyn[name], ctx.out_filename, name, start=start, count=cnt)
            end

            if iteration == 1 # Only write static outputs once
                for name ∈ ctx.out_vars_stat
                    if name == ctx.stat_var_mapping["b"] # Water height correction
                        itbuf.out_grids_stat[name] .-= ctx.water_height
                    end
                    ncwrite(itbuf.out_grids_stat[name], ctx.out_filename, name, start=[1, 1], count=[-1, -1])
                end
            end

            println("Done.")
        end # for iteration
    end

    function _bins_equal_load(ctx)
        # Constant load factor per column
        return ones(Float64, ctx.samples[X])
    end

    function bin(ctx)
        println("Binning $(ctx.n_simplices) $(simplex_name(ctx)) into $(ctx.n_threads)+$(ctx.n_threads - 1)+1 buckets...")

        # Assign a load factor to each x-column. Currently just 1.0 for each column, but dynamically calculating factors based on the simplices in each
        # column might also be an option.
        l_x_simps = _bins_equal_load(ctx)

        # Now that each column has a load factor, size each bin such that the load is distributed evenly.
        l_bin_start_idxs = Array{Int, 1}(undef, ctx.n_threads)
        n_opt_simps_per_p_bin = sum(l_x_simps) / ctx.n_threads
        l_bin_start_idxs[1] = 1

        # Start from the left bin and add columns to it until its rounded load factor is closest to the optimal value (1/n_threads*total_load)
        # Then continue with next bin
        next_x_idx = 1
        for bin_id ∈ 1:(ctx.n_threads-1)
            current_simp_count = 0.
            new_simp_count = 0.
            while (n_opt_simps_per_p_bin - current_simp_count > new_simp_count - n_opt_simps_per_p_bin) && next_x_idx <= ctx.samples[X]
                current_simp_count = new_simp_count
                new_simp_count = current_simp_count + l_x_simps[next_x_idx]
                next_x_idx += 1
            end

            l_bin_start_idxs[bin_id + 1] = min(next_x_idx, ctx.samples[X])
        end

        #============================================#
        # Perform binning
        #============================================#

        # For each tetrahedron, assign a bin with an ID ∈ 1:2*n_threads based on its location
        l_bin_ids = Array{UInt16, 1}(undef, ctx.n_simplices)

        Threads.@threads for thread_id ∈ 1:ctx.n_threads
            thread_range_min = floor(Int, (thread_id - 1) * ctx.simplices_per_thread) + 1
            thread_range_max = thread_id == ctx.n_threads ? ctx.n_simplices : floor(Int, thread_id * ctx.simplices_per_thread)

            m43_tet_points = Array{Float64, 2}(undef, (ctx.n_simplex_points, 3))

            for tet_id ∈ thread_range_min:thread_range_max
                tetrahedron = (@view ctx.simplices[:, tet_id])
                for i ∈ 1:ctx.n_simplex_points
                    m43_tet_points[i, :] = ctx.points[:, tetrahedron[i]]
                end

                # Simplex is not in the z-domain (only possible if z_range ∈ [z_floor, z_surface]), ignore it.
                if ctx.z_range != z_all
                    is_bath = is_bathy(m43_tet_points, ctx.water_height)

                    if (ctx.z_range == z_floor) != is_bath
                        l_bin_ids[tet_id] = ctx.n_threads * 2 + 1 # Unused bin, will be ignored later
                        continue
                    end
                else # This would be were tetrahedra in the volume mesh that belong to the ground are filtered out. For now: cutoff below -2000m
                    # TODO: correct z filter for tetrahedra
                    coord_z_min = minimum(m43_tet_points[:, 3])
                    coord_z_max = maximum(m43_tet_points[:, 3])
                    if coord_z_max <= 0.00001
                        l_bin_ids[tet_id] = ctx.n_threads * 2 + 1 # Unused bin, will be ignored later
                        continue
                    end
                end

                coord_x_min = minimum(m43_tet_points[:, 1])
                coord_x_max = maximum(m43_tet_points[:, 1])

                idx_x_min = floor(Int,(coord_x_min - ctx.domain[X][MIN]) / ctx.sampling_rate[X]) + 1
                idx_x_max =  ceil(Int,(coord_x_max - ctx.domain[X][MIN]) / ctx.sampling_rate[X])

                bin_id_l = findlast(bin_start_idx -> bin_start_idx ≤ idx_x_min, l_bin_start_idxs)
                bin_id_r = findlast(bin_start_idx -> idx_x_max ≥ bin_start_idx, l_bin_start_idxs)
                if bin_id_r > ctx.n_threads; bin_id_r = ctx.n_threads; end # The right border of the last bin is inclusive while the others are not.

                # n_bins := 2*n_threads
                # bin i ∈ 1:n_bins contains the tets completely inside the column of bin i
                # bin j ∈ n_bins+1:2*n_bins-1 contains tets that are overlapping across bin i and i+1
                # bin k = 2*n_bins contains tets overlapping over more than 2 adjascent columns
                if bin_id_l == bin_id_r
                    l_bin_ids[tet_id] = bin_id_l
                elseif bin_id_l == bin_id_r - 1
                    l_bin_ids[tet_id] = ctx.n_threads + bin_id_l
                else
                    l_bin_ids[tet_id] = 2 * ctx.n_threads
                end
            end # for tet_id
        end # for thread_id

        #============================================#
        # Calculate number of simps in each bin,
        # print statistics
        #============================================#

        # The number of tetrahedra in each respective bin (only the first 2*n_threads bins have been filled)
        l_bin_counts = zeros(Int, 2 * ctx.n_threads + 1)

        for bin_id ∈ l_bin_ids
            l_bin_counts[bin_id] += 1
        end

        for bin_id ∈ 1:length(l_bin_counts) - 1
            println(l_bin_counts[bin_id], " ", simplex_name(ctx), " in bin ", bin_id)
        end

        if ctx.z_range != z_all
            println(l_bin_counts[length(l_bin_counts)], " ", simplex_name(ctx), " are ignored due to z_range filter.")
        end

        println("Done.")
        return (l_bin_ids, l_bin_counts)
    end

    function iterate_simps!(bin_id, ctx:: RasterizationContext, itbuf:: IterationBuffers, t_start, t_stop; print_progress=false)
        n_times = t_stop - t_start + 1

        #============================================#
        # Pre-allocate memory used throughout
        # rasterization by this thread
        #============================================#

        m43_simp_points     = Array{Float64, 2}(undef, (ctx.n_simplex_points, 3))
        m32_tet_aabb        = Array{Float64, 2}(undef, (3, 2))
        l_face_points       = Array{SubArray, 1}(undef, ctx.n_simplex_points)
        l_tet_faces         = Array{HesseNormalForm, 1}(undef, ctx.n_simplex_points)

        v_p                 = Array{Float64, 1}(undef, 3)
        v_n                 = Array{Float64, 1}(undef, 3)
        simp_sample_counts  = Array{UInt16, 2}(undef, (64, 64))

        d_start             = now()
        d_last_printed      = d_start
        print_interval      = Second(2)
        n_bin_rasterized    = 0
        n_bin_total         = itbuf.bin_counts[bin_id]

        for simp_id ∈ 1:ctx.n_simplices
            # Only process tetrahedra in thread's own bin
            if itbuf.bin_ids[simp_id] != bin_id
                continue
            end

            if print_progress
                d_now = now()

                if (d_now - d_last_printed) >= print_interval && n_bin_rasterized > 0
                    d_last_printed = d_now

                    etr  = (d_now-d_start) * (n_bin_total - n_bin_rasterized) ÷ (n_bin_rasterized)
                    hh   = floor(etr, Dates.Hour)
                    etr -= floor(hh,  Dates.Millisecond)
                    mm   = floor(etr, Dates.Minute)
                    etr -= floor(mm,  Dates.Millisecond)
                    ss   = floor(etr, Dates.Second)

                    @printf("Working on simplex %7d of %7d in timesteps %d-%d (%02.2f%% done). ETR: %02d:%02d:%02d\n",
                            n_bin_rasterized, n_bin_total,
                            t_start - 1,           # Tools like ParaView work 0-indexed. Avoid confusion by outputting 0-indexed here.
                            t_stop  - 1,           # Tools like ParaView work 0-indexed. Avoid confusion by outputting 0-indexed here.
                            n_bin_rasterized / n_bin_total * 100,
                            hh.value, mm.value, ss.value)
                end
            end # if print_progress

            # Read simplex points. m43 is only filled with 3x3 values for triangles
            simplex_point_idxs = (@view ctx.simplices[:, simp_id])
            for i ∈ 1:ctx.n_simplex_points
                m43_simp_points[i, X:Z] = ctx.points[X:Z, simplex_point_idxs[i]]
            end

            for dim ∈ X:Z
                m32_tet_aabb[dim, MIN] = minimum(m43_simp_points[:, dim])
                m32_tet_aabb[dim, MAX] = maximum(m43_simp_points[:, dim])
            end

            #============================================#
            # A note on global & local indices
            #============================================#
            # In the following, idx_g_ will refer to an
            # index of the global sampling domain while
            # idx_l_ refers to an index in the simp-local
            # domain.
            #============================================#

            idx_g_min = Tuple(
                        (floor(Int32, (m32_tet_aabb[dim,MIN] - ctx.domain[dim][MIN]) / ctx.sampling_rate[dim]) + 1) for dim ∈ X:Z)

            idx_g_max = ( ceil(Int32, (m32_tet_aabb[  X,MAX] - ctx.domain[  X][MIN]) / ctx.sampling_rate[  X]),
                          ceil(Int32, (m32_tet_aabb[  Y,MAX] - ctx.domain[  Y][MIN]) / ctx.sampling_rate[  Y]),
                         floor(Int32, (m32_tet_aabb[  Z,MAX] - ctx.domain[  Z][MIN]) / ctx.sampling_rate[  Z]) + 1)

            n_current_cells = idx_g_max .- idx_g_min .+ 1

            # If the current simp has a bounding box larger than anticipated, enlarge the sample counts array accordingly
            # By allocating double the size of the bounding box, future simps that are just a bit larger than the current one
            # will not trigger a resize operation.
            if size(simp_sample_counts, X) < n_current_cells[X] || size(simp_sample_counts, Y) < n_current_cells[Y]
                simp_sample_counts = Array{UInt16, 2}(undef, (n_current_cells[X] * 2, n_current_cells[Y] * 2))
            end

            # For the bounding box of the current tet, reset all sample counts to 0
            simp_sample_counts[1:n_current_cells[X], 1:n_current_cells[Y]] .= 0

            for i ∈ 1:ctx.n_simplex_points
                l_face_points[i] = @view m43_simp_points[i,:]
            end

            for i ∈ 1:ctx.n_simplex_points
                # Exclude point i from the tetrahedron and consider the plane defined by the other 3 points
                # If and only if (x,y,z) lies on the "positive" side (in normal direction) of all 4 planes,
                # it is inside the tetrahedron.
                # Analogous for triangles but the fourth point is placed straight above the third one which
                # allows the usage of the same point-in-simplex code afterwards.
                v_p_excl = l_face_points[i]
                v_p1 =     l_face_points[( i    % ctx.n_simplex_points) + 1]
                v_p2 =     l_face_points[((i+1) % ctx.n_simplex_points) + 1]

                v_p3 =
                    if ctx.n_simplex_points == 4; l_face_points[((i+2) % 4) + 1]
                    elseif ctx.n_simplex_points == 3; [v_p1[1], v_p1[2], v_p1[3] + 1.]
                    else throw(DimensionMismatch("Only simplices with 3 or 4 points (read: triangles / tetrahedra) are supported."))
                    end

                # Calculate Hesse normal form of the plane defined by p1, p2, p3
                # Normal vector
                cross3!(v_p2-v_p1, v_p3-v_p1, v_n)
                ctx.n_simplex_points == 3 && @assert (v_n[3] == 0.) "Normal vector of 2D line has non-zero z component!"
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
                # For each tetrahedron face, filter out cells that are not in the tetrahedron
                for i ∈ 1:ctx.n_simplex_points
                    face = l_tet_faces[i]
                    # p_excl is ALWAYS on the "positive" side of the plane.
                    # So, if p_excl and p lie on the same side, p is inside the tetrahedron
                    dist_p = dot3(v_p, face.n0) + face.d

                    # If signs of distance vectors do not match, p and p_excl are on different sides of the plane (reject)
                    # As stated above, the normal vector n is chosen so that p_excl is on the POSITIVE side of the plane
                    # Therefore, we only need to check if dist_p is negative to reject p.
                    if dist_p < -EDGE_TOL
                        return false
                    end

                    # The point lies approx. ON the face, reject in the tet on the lexicographically
                    # more negative side of the face. This ensures that for two neighboring tetrahedra, the point is
                    # only accepted for one of them.
                    # Def.: lexicographical order:
                    #       a ≤ b ⟺ (a.x ≤ b.x) ∨ (a.x = b.x ∧ a.y ≤ b.y) ∨ (a.x = b.x ∧ a.y = b.y ∧ a.z ≤ b.z)
                    # This is a total ordering on all vectors in R^3 and will therefore always be larger for the "positive"
                    # normal vector in contrast to its negated counterpart.
                    if (abs(dist_p) <= EDGE_TOL)
                        is_positive_side =
                            (face.n0[X]  > 0.) ||
                            (face.n0[X] == 0. && face.n0[Y]  > 0.) ||
                            (face.n0[X] == 0. && face.n0[Y] == 0. && face.n0[Z] > 0.)

                        # Reject if in tetrahedron on face's lex. negative side
                        if !is_positive_side
                            return false
                        end
                    end
                end # for i

                return true
            end

            #============================================#
            # Sample each cuboid that is overlapping the
            # current simplex's AABB.
            # If its CENTER POINT lies within the simp,
            # increase the simp's sample count for the
            # cuboid's (x, y) index by 1.
            #============================================#

            for idx_g_x ∈ idx_g_min[X]:idx_g_max[X]
                v_p[1] = (idx_g_x-1) * ctx.sampling_rate[X] + ctx.domain[X][MIN]

                for idx_g_y ∈ idx_g_min[Y]:idx_g_max[Y]
                    v_p[2] = (idx_g_y-1) * ctx.sampling_rate[Y] + ctx.domain[Y][MIN]

                    if ctx.n_simplex_points == 4 # processing tets, increment counter for intersected cells

                        for idx_g_z ∈ idx_g_min[Z]:idx_g_max[Z]
                            v_p[3] = (idx_g_z-1) * ctx.sampling_rate[Z] + ctx.domain[Z][MIN]

                            if is_in_simplex() # accept cell and increment its 2D cell counter for the current tetrahedron by 1
                                idx_l_x = idx_g_x - idx_g_min[X] + 1
                                idx_l_y = idx_g_y - idx_g_min[Y] + 1
                                itbuf.out_sample_counts[idx_l_x, idx_l_y] += 1
                            end
                        end # for z
                    else # processing triangles, immediately rasterize vars
                        if is_in_simplex()
                            idx_l_x = idx_g_x - idx_g_min[X] + 1
                            idx_l_y = idx_g_y - idx_g_min[Y] + 1

                            # Only calculate static vars in the first timestep (read: once)
                            if ctx.n_out_vars_stat != 0 && t_start == ctx.t_begin
                                v_λ = bary_coords_2d(v_p, m43_simp_points)

                                # Calculate bathymetry height from triangle's z-coords and interpolate between them
                                b = sum(m43_simp_points[1:3,Z] .* v_λ)
                                itbuf.out_grids_stat[ctx.stat_var_mapping["b"]][idx_g_x, idx_g_y] = b
                            end

                            if ctx.tanioka
                                tanioka_offset = 0.

                                # ∂b/∂x
                                v_λ = bary_coords_2d(v_p .+ (1., 0., 0.), m43_simp_points)
                                b = sum(m43_simp_points[1:3,Z] .* v_λ)
                                ∂b∂x = b - itbuf.out_grids_stat[ctx.stat_var_mapping["b"]][idx_g_x, idx_g_y]

                                # ∂b/∂y
                                v_λ = bary_coords_2d(v_p .+ (0., 1., 0.), m43_simp_points)
                                b = sum(m43_simp_points[1:3,Z] .* v_λ)
                                ∂b∂y = b - itbuf.out_grids_stat[ctx.stat_var_mapping["b"]][idx_g_x, idx_g_y]

                                for t ∈ 1:n_times
                                    # TODO: Do not hardcode variable indices for U,V,W!!
                                    tanioka_offset = (∂b∂x * itbuf.prefetched_vars["U"][simp_id, t] 
                                                    +  ∂b∂y * itbuf.prefetched_vars["V"][simp_id, t])
                                    itbuf.out_grids_dyn[ctx.dyn_var_mapping["W"]][idx_g_x, idx_g_y, t] = itbuf.prefetched_vars["W"][simp_id, t] - tanioka_offset
                                end
                            else
                                for in_name ∈ keys(ctx.dyn_var_mapping)
                                    out_name = ctx.dyn_var_mapping[in_name]
                                    for t ∈ 1:n_times
                                        itbuf.out_grids_dyn[out_name][idx_g_x, idx_g_y, t] = itbuf.prefetched_vars[in_name][simp_id, t]
                                    end
                                end
                            end

                            itbuf.out_sample_counts[idx_l_x, idx_l_y] = 1 # 2D surface is only sampled once in z-direction
                        end # if is_in_simplex()
                    end # if ctx.n_simplex_points == 4
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

            if ctx.z_range == z_all # Build running average of variables in volume grid
                for idx_l_y ∈ 1:n_current_cells_y, idx_l_x ∈ 1:n_current_cells_x
                    num_samples = itbuf.out_sample_counts[idx_l_x, idx_l_y]

                    if num_samples > 0
                        idx_g_x = idx_l_x + idx_g_min[X] - 1
                        idx_g_y = idx_l_y + idx_g_min[Y] - 1
                        total_samples = (itbuf.out_sample_counts[idx_g_x, idx_g_y] += num_samples)

                        for in_name ∈ keys(ctx.dyn_var_mapping)
                            out_name = ctx.dyn_var_mapping[in_name]

                            for t ∈ 1:n_times
                                avg = itbuf.out_grids_dyn[out_name][idx_g_x, idx_g_y, t]
                                itbuf.out_grids_dyn[out_name][idx_g_x, idx_g_y, t] = 
                                    avg + (itbuf.prefetched_vars[in_name][simp_id, t] - avg) * (num_samples / total_samples)
                            end
                        end
                    end
                end
            end # if z_range

            n_bin_rasterized += 1
        end # for tet_id
    end

    """
    Checks if the given triangle belongs to the bathymetry mesh.
        This is done as follows:
        If the triangle is not parallel to the xy-plane (different z-coords for at least one of the points), it belongs to bathymetry.
        If the triangle is parallel to the xy-plane but not at water height, then it also belongs to bathymetry.
        Otherwise, it does not belong to bathymetry.

    # Arguments
    - `triangle_points`:    The 2D array of triangle point coords (each column consists of x, y, and z coord)
    - `water_height`:       The water height as a number
    """
    @inline function is_bathy(triangle_points, water_height = 0.)
        min = minimum(triangle_points[:, 3])
        if abs(min - water_height) > 1e-8; return true; end # Check height condition

        max = maximum(triangle_points[:, 3])
        if abs(max - min) > 1e-8; return true; end # Check flatness condition

        return false
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

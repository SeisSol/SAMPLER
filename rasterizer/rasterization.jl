#=
rasterization.jl:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

module Rasterization
    using Dates: length
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

    """
    Contains all static information necessary for rasterization, as well as some computed values for convenience.
    """
    struct RasterizationContext
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

    """
    Contains buffers that are rewritten during rasterization, as well as binning info.
    """
    struct IterationBuffers
        bin_ids
        bin_counts
        out_grids_dyn               :: Dict{String, Array{Float64, 3}}
        out_grids_stat              :: Dict{String, Array{Float64, 2}}
        out_sample_counts           :: Array{UInt16, 2}
        prefetched_vars             :: Dict{String, Array{Float64, 2}}
    end

    """
    Returns the name of the simplex type currently being rasterized (triangle or tetrahedron), either in plural (default) or singular form.
    """
    @inline function simplex_name(ctx, plural=true)
        return ctx.n_dims == 3 ? (plural ? "simplices" : "simplex") : (plural ? "triangles" : "triangle") 
    end

    function RasterizationContext(simplices, points, sampling_rate, mem_limit, dyn_var_mapping, stat_var_mapping, z_range, 
                                  tanioka, water_height, t_begin, t_end, times, out_filename, custom_domain; max_threads=typemax(Int))

        custom_domain                   = (custom_domain..., (-Inf, Inf)) # Add z-range
        domain                          = Tuple( (max(minimum(points[i, :]), custom_domain[i][MIN]), 
                                                  min(maximum(points[i, :]), custom_domain[i][MAX]))                for i in X:Z)

        samples                         = Tuple( ceil(UInt, (domain[i][MAX] - domain[i][MIN]) / sampling_rate[i])   for i in X:Z)

        #============================================#
        # Calculate problem size & hardware limitations
        #============================================#

        n_threads                       = min(max_threads, nthreads())

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

        # static vars are currently only bathymetry and thus extracted from grid geometry. Thus, no in_vars_stat.
        in_vars_dyn                     = collect(keys(dyn_var_mapping))

        # Tanioka's method also needs horizontal displacements, add them if not present yet. 
        if tanioka
            for var_name ∈ ("U", "V")
                var_name ∈ in_vars_dyn || push!(in_vars_dyn, var_name)
            end
        end

        n_in_vars                       = length(in_vars_dyn)

        out_vars_stat                   = z_range != z_all ? collect(values(stat_var_mapping)) : []
        out_vars_dyn                    = collect(values(dyn_var_mapping))

        n_out_vars_dyn                  = length(dyn_var_mapping)
        n_out_vars_stat                 = length(out_vars_stat)

        n_timesteps                     = t_end - t_begin + 1
        n_timesteps_per_iteration       = ((mem_limit - b_mem_points - b_mem_simps - b_mem_cnt - b_mem_misc - n_out_vars_stat * b_mem_per_out_var_and_timestep)
                                        ÷  (n_in_vars * b_mem_per_in_var_and_timestep + n_out_vars_dyn * b_mem_per_out_var_and_timestep))
        n_timesteps_per_iteration       = max(1, n_timesteps_per_iteration) # Process at least 1 timestep

        n_iterations                    = ceil(Int, n_timesteps / n_timesteps_per_iteration)

        return RasterizationContext(simplices, points, 
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
    - `xdmf`:                   The XDMF file handle of the file to be rasterized
    - `dyn_var_mapping`:        Mapping of input to output names for dynamic variables in the current rasterization step.
    - `stat_var_mapping`:       Analogous, but for static variables.
    - `times`:                  An array of timestamps corresponding to the timesteps
    - `sampling_rate`:          A tuple (Δx, Δy, Δz) defining the regular grid cell size.
    - `out_filename`:           The path/name of the output NetCDF file
    - `mem_limit`:              The soft memory (in bytes) the software is allowed to utilize (soft limit!)
    - `z_range`:                A filter for discarding simplices in a certain mesh region. `z_floor` only keeps bathymetry, `z_surface` only keeps the sea surface, 
                                `z_all` keeps all simplices.
    - `t_begin`:                The first timestep to be rasterized
    - `t_end`:                  The last timestep to be rasterized
    - `create_nc_dyn_vars`:     If not empty: create output file containing all value sides of the var mapping as variables. Existing file is overwritten.
                                Else: use existing file. If given, should contain all variable mappings of all rasterization passes.
    - `create_nc_stat_vars`:    Analogous to above.
    - `water_height`:           The offset to translate the water level to be at `z=0`.
    - `tanioka`:                Whether to apply Tanioka's method to the bathymetry. Requires U, V, W displacements to be present in the input XDMF file.
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
        tanioka             :: Bool    = false,
        domain              :: Main.Args.DomainSize = ((-Inf, Inf),(-Inf, Inf)))

        simplices, points = grid_of(xdmf)

        ctx = RasterizationContext(simplices, points, sampling_rate, mem_limit, dyn_var_mapping, stat_var_mapping, z_range, tanioka, 
                                   water_height, t_begin, t_end, times, out_filename, domain)

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

        iterate(ctx, itbuf, xdmf)

    end # function rasterize

    function bin(ctx)
        println("Binning $(ctx.n_simplices) $(simplex_name(ctx)) into $(ctx.n_threads)+$(ctx.n_threads - 1)+1 buckets...")

        n_excluded_domain = 0
        n_excluded_z_range = 0

        l_bin_start_idxs = [1 + round(Int, ctx.samples[X] / ctx.n_threads * bucket_id) for bucket_id ∈ 0:(ctx.n_threads-1)]

        #============================================#
        # Perform binning
        #============================================#

        # For each tetrahedron, assign a bin with an ID ∈ 1:2*n_threads based on its location
        l_bin_ids = Array{UInt16, 1}(undef, ctx.n_simplices)

        Threads.@threads for thread_id ∈ 1:ctx.n_threads
            thread_range_min = floor(Int, (thread_id - 1) * ctx.simplices_per_thread) + 1
            thread_range_max = thread_id == ctx.n_threads ? ctx.n_simplices : floor(Int, thread_id * ctx.simplices_per_thread)

            m3n_simp_points = Array{Float64, 2}(undef, (3, ctx.n_simplex_points))

            for simp_id ∈ thread_range_min:thread_range_max
                tetrahedron = (@view ctx.simplices[:, simp_id])
                for i ∈ 1:ctx.n_simplex_points
                    m3n_simp_points[:, i] = ctx.points[:, tetrahedron[i]]
                end

                coord_x_min = minimum(m3n_simp_points[X, :])
                coord_x_max = maximum(m3n_simp_points[X, :])
                coord_y_min = minimum(m3n_simp_points[Y, :])
                coord_y_max = maximum(m3n_simp_points[Y, :])

                if coord_x_max <= ctx.domain[X][MIN] || coord_x_min >= ctx.domain[X][MAX] || coord_y_max <= ctx.domain[Y][MIN] || coord_y_min >= ctx.domain[Y][MAX]
                    l_bin_ids[simp_id] = ctx.n_threads * 2 + 1
                    n_excluded_domain += 1
                    continue
                end

                # If simplex is not in the z-domain (only possible if z_range ∈ [z_floor, z_surface]), ignore it.
                if ctx.z_range != z_all
                    is_bath = is_bathy(m3n_simp_points, ctx.water_height)

                    if (ctx.z_range == z_floor) != is_bath
                        l_bin_ids[simp_id] = ctx.n_threads * 2 + 1 # Unused bin, will be ignored later
                        n_excluded_z_range += 1
                        continue
                    end
                else # This would be were tetrahedra in the volume mesh that belong to the ground are filtered out. For now: cutoff below -2000m
                    # TODO: correct z filter for tetrahedra
                    coord_z_min = minimum(m3n_simp_points[Z, :])
                    coord_z_max = maximum(m3n_simp_points[Z, :])
                    if coord_z_max <= 0.00001
                        l_bin_ids[simp_id] = ctx.n_threads * 2 + 1 # Unused bin, will be ignored later
                        n_excluded_z_range += 1
                        continue
                    end
                end

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
                    l_bin_ids[simp_id] = bin_id_l
                elseif bin_id_l == bin_id_r - 1
                    l_bin_ids[simp_id] = ctx.n_threads + bin_id_l
                else
                    l_bin_ids[simp_id] = 2 * ctx.n_threads
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

        max_chars_per_bin = floor(Int, max(0, log10(maximum(l_bin_counts)))) + 1
        chars_per_line = 120
        bins_per_line = max(1, (chars_per_line ÷ (max_chars_per_bin + 1)) ÷ 8 * 8)

        println("Bin Counts:")
        for bin_id ∈ 1:length(l_bin_counts) - 1
            if bin_id != 1 && (bin_id - 1) % bins_per_line == 0
                println()
            end
            bin_chars = floor(Int, max(0, log10(l_bin_counts[bin_id]))) + 1
            print(' '^(max_chars_per_bin-bin_chars+1), l_bin_counts[bin_id])
        end
        println()

        if n_excluded_domain != 0
            println(n_excluded_domain, " ", simplex_name(ctx), " are ignored due to domain filter.")
        end

        if n_excluded_z_range != 0
            println(n_excluded_z_range, " ", simplex_name(ctx), " are ignored due to z_range filter.")
        end

        println("Done.")
        return (l_bin_ids, l_bin_counts)
    end

    function iterate(ctx, itbuf, xdmf)
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

            var_bufs :: Dict{AbstractString, Array{AbstractArray, 1}} = Main.XDMF.data_of(xdmf, t_start, t_stop, ctx.in_vars_dyn)

            for var_name ∈ ctx.in_vars_dyn, t ∈ 1:n_times
                copyto!(itbuf.prefetched_vars[var_name], 1 + (t-1) * ctx.n_simplices, var_bufs[var_name][t], 1, ctx.n_simplices)
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

    function iterate_simps!(bin_id, ctx:: RasterizationContext, itbuf:: IterationBuffers, t_start, t_stop; print_progress=false)
        n_times = t_stop - t_start + 1

        #============================================#
        # Pre-allocate memory used throughout
        # rasterization by this thread
        #============================================#

        m43_simp_points     = Array{Float64, 2}(undef, (3, ctx.n_simplex_points))
        m32_tet_aabb        = Array{Float64, 2}(undef, (3, 2))
        l_simp_points       = Array{SubArray, 1}(undef, ctx.n_simplex_points)
        l_simp_hnfs         = Array{HesseNormalForm, 1}(undef, ctx.n_simplex_points)

        v_p                 = Array{Float64, 1}(undef, 3)
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
            for i ∈ 1:ctx.n_simplex_points
                m43_simp_points[:, i] = ctx.points[:, ctx.simplices[i, simp_id]]
            end

            for dim ∈ X:Z
                m32_tet_aabb[dim, MIN] = minimum(m43_simp_points[dim, :])
                m32_tet_aabb[dim, MAX] = maximum(m43_simp_points[dim, :])
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

            idx_g_min = Tuple(max(1               , idx_g_min[dim]) for dim ∈ X:Z)
            idx_g_max = Tuple(min(ctx.samples[dim], idx_g_max[dim]) for dim ∈ X:Z)

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
                l_simp_points[i] = @view m43_simp_points[:,i]
            end

            calculate_hnfs!(ctx, l_simp_points, l_simp_hnfs)

            #============================================#
            # Sample each cuboid that is overlapping the
            # current simplex's AABB.
            # If its CENTER POINT lies within the simp,
            # increase the simp's sample count for the
            # cuboid's (x, y) index by 1.
            #============================================#

            for idx_g_x ∈ idx_g_min[X]:idx_g_max[X]
                v_p[X] = (idx_g_x-1) * ctx.sampling_rate[X] + ctx.domain[X][MIN]

                for idx_g_y ∈ idx_g_min[Y]:idx_g_max[Y]
                    v_p[Y] = (idx_g_y-1) * ctx.sampling_rate[Y] + ctx.domain[Y][MIN]

                    if ctx.n_simplex_points == 4 # processing tets, increment counter for intersected cells

                        for idx_g_z ∈ idx_g_min[Z]:idx_g_max[Z]
                            v_p[Z] = (idx_g_z-1) * ctx.sampling_rate[Z] + ctx.domain[Z][MIN]

                            if is_in_simplex(ctx, l_simp_hnfs, v_p) # accept cell and increment its 2D cell counter for the current tetrahedron by 1
                                idx_l_x = idx_g_x - idx_g_min[X] + 1
                                idx_l_y = idx_g_y - idx_g_min[Y] + 1
                                itbuf.out_sample_counts[idx_l_x, idx_l_y] += 1
                            end
                        end # for z
                    else # processing triangles, immediately rasterize vars
                        if is_in_simplex(ctx, l_simp_hnfs, v_p)
                            idx_l_x = idx_g_x - idx_g_min[X] + 1
                            idx_l_y = idx_g_y - idx_g_min[Y] + 1

                            # Only calculate static vars in the first timestep (read: once)
                            if ctx.n_out_vars_stat != 0 && t_start == ctx.t_begin
                                v_λ = bary_coords_2d(v_p, m43_simp_points)

                                # Calculate bathymetry height from triangle's z-coords and interpolate between them
                                b = sum(m43_simp_points[Z,1:3] .* v_λ)
                                itbuf.out_grids_stat[ctx.stat_var_mapping["b"]][idx_g_x, idx_g_y] = b
                            end

                            if ctx.tanioka
                                (∂b∂x, ∂b∂y) = get_tanioka_db(ctx, itbuf, v_p, m43_simp_points, idx_g_x, idx_g_y)
                            else
                                (∂b∂x, ∂b∂y) = (0., 0.)
                            end

                            for in_name ∈ keys(ctx.dyn_var_mapping)
                                out_name = ctx.dyn_var_mapping[in_name]
                                for t ∈ 1:n_times
                                    if in_name == "W" && ctx.tanioka
                                        tanioka_offset = (∂b∂x * itbuf.prefetched_vars["U"][simp_id, t] +  ∂b∂y * itbuf.prefetched_vars["V"][simp_id, t])
                                    else
                                        tanioka_offset = 0.
                                    end
                                    
                                    itbuf.out_grids_dyn[out_name][idx_g_x, idx_g_y, t] = itbuf.prefetched_vars[in_name][simp_id, t] - tanioka_offset
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
        Calculates all 3 / 4 Hesse Normal Forms for the given simplex points.

    # Arguments
    - `ctx` The rasterization context
    - `l_simp_points` a list of ctx.n_simplex_points column vectors.
    - `l_simp_hnfs` a HesseNormalForm-array of size ctx.n_simplex_points. Will be populated with the calculated HNFs.
    """
    @inline function calculate_hnfs!(ctx, l_simp_points, l_simp_hnfs)
        v_n = Array{Float64, 1}(undef, 3)

        for i ∈ 1:ctx.n_simplex_points
            # Exclude point i from the tetrahedron and consider the plane defined by the other 3 points
            # If and only if (x,y,z) lies on the "positive" side (in normal direction) of all 4 planes,
            # it is inside the tetrahedron.
            # Analogous for triangles but the fourth point is placed straight above the third one which
            # allows the usage of the same point-in-simplex code afterwards.
            v_p_excl = l_simp_points[i]
            v_p1     = l_simp_points[( i    % ctx.n_simplex_points) + 1]
            v_p2     = l_simp_points[((i+1) % ctx.n_simplex_points) + 1]

            v_p3 =
                if ctx.n_simplex_points == 4; l_simp_points[((i+2) % 4) + 1]
                elseif ctx.n_simplex_points == 3; [v_p1[1], v_p1[2], v_p1[3] + 1.]
                else throw(DimensionMismatch("Only simplices with 3 or 4 points (read: triangles / tetrahedra) are supported."))
                end

            # Calculate Hesse normal form of the plane defined by p1, p2, p3
            # Normal vector
            cross3!(v_p2-v_p1, v_p3-v_p1, v_n)
            ctx.n_simplex_points == 3 && @assert (v_n[3] == 0.) "Normal vector of 2D line has non-zero z component!"
            normalize!(v_n)

            # n1x1 + n2x2 + n3x3 - d = 0; solve for d
            d = dot3(v_p1, v_n)
            dist_p_excl = dot3(v_p_excl, v_n) - d

            # Ensure that p_excl is considered to be on the "positive" side of the plane
            if dist_p_excl < 0
                v_n = -v_n
                d = -d
            end

            l_simp_hnfs[i] = HesseNormalForm((v_n[1], v_n[2], v_n[3]), d)
        end # for i
    end

    @inline function is_in_simplex(ctx, l_tet_faces, point)
        # For each tetrahedron face, filter out cells that are not in the tetrahedron
        for i ∈ 1:ctx.n_simplex_points
            face = l_tet_faces[i]
            # p_excl is ALWAYS on the "positive" side of the plane.
            # So, if p_excl and p lie on the same side, p is inside the tetrahedron
            dist_p = dot3(point, face.n0) - face.d

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
        return any(abs.(triangle_points[Z,:] .- water_height) .> 1e-8)
    end

    @inline function get_tanioka_db(ctx, itbuf, v_p, m43_simp_points, idx_g_x, idx_g_y)
        # ∂b/∂x
        v_λ = bary_coords_2d(v_p .+ (1., 0., 0.), m43_simp_points)
        b = sum(m43_simp_points[Z,1:3] .* v_λ)
        ∂b∂x = b - itbuf.out_grids_stat[ctx.stat_var_mapping["b"]][idx_g_x, idx_g_y]

        # ∂b/∂y
        v_λ = bary_coords_2d(v_p .+ (0., 1., 0.), m43_simp_points)
        b = sum(m43_simp_points[Z,1:3] .* v_λ)
        ∂b∂y = b - itbuf.out_grids_stat[ctx.stat_var_mapping["b"]][idx_g_x, idx_g_y]

        return (∂b∂x, ∂b∂y)
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
        @assert size(ps, 1) >= 2
        @assert size(ps, 2) >= 3

        T11 = ps[X, 1] - ps[X, 3]
        T12 = ps[X, 2] - ps[X, 3]
        T21 = ps[Y, 1] - ps[Y, 3]
        T22 = ps[Y, 2] - ps[Y, 3]

        det = T11 * T22 - T12 * T21

        λ1 = (T22*(p[X]-ps[X,3]) - T12*(p[Y]-ps[Y,3])) / det
        λ2 = (T11*(p[Y]-ps[Y,3]) - T21*(p[X]-ps[X,3])) / det
        λ3 = 1. - λ1 - λ2

        return (λ1, λ2, λ3)
    end
end # module Rasterization

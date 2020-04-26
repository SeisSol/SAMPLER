module Rasterization
    using Base.Threads
    using LinearAlgebra
    using Dates
    using Printf
    using NetCDF
    using Profile
    using Main.Util

    struct TetrahedronFace
        n               :: NTuple{3, Float64}
        d               :: Float64
        dist_p_excl     :: Float64
    end

    mutable struct ThreadContext
        tet_points      :: Array{Float64, 2}
        tet_aabb        :: Array{Float64, 2}
        tet_dict        :: Array{Int32, 2}
        p               :: Array{Float64, 1}
        face_points     :: Array{SubArray, 1}
        n               :: Array{Float64, 1}
        tet_faces       :: Array{TetrahedronFace, 1}
    end

    const EDGE_TOL = .0

    #=
    tetrahedra, points_3d, XDMF.data_of(args["input-file-3d"], var_names_3d...), var_names_3d, 
    triangles,  points_2d, XDMF.data_of(args["input-file-2d"], (var_names_floor+var_names_surface)...), 
    var_names_floor, var_names_surface, times_3d, sampling_rate, 
    args["output-file"], args["memory-limit"]
    =#

    function rasterize(
        tetrahedra          :: AbstractArray{INDEX_TYPE, 2}, 
        points_3d           :: AbstractArray{Float64, 2},
        vars_3d             :: AbstractArray{A where A <: AbstractArray, 2},
        var_names_3d        :: AbstractArray,
        triangles           :: AbstractArray{INDEX_TYPE, 2}, 
        points_2d           :: AbstractArray{Float64, 2},
        vars_2d             :: AbstractArray{A where A <: AbstractArray, 2}, 
        var_names_floor     :: AbstractArray,
        var_names_surface   :: AbstractArray,
        times               :: AbstractArray{Float64, 1}, 
        sampling_rate       :: NTuple{3, Float64}, 
        out_filename        :: AbstractString,
        mem_limit           :: Int64) where INDEX_TYPE <: Integer

        #============================================#
        # Code style & conventions
        #============================================#
        #
        # Some variables are prefixed with <prefix>_
        # to indicate their type or unit:
        #   i_  : Interval / Pair / 2-Tuple
        #   n_  : Count / Number of Objects
        #   b_  : Size in Bytes
        #   l_  : 1D-Array / List
        #   d_  : Date & Time
        #   mMN_: M×N Matrix
        #   vN_ : N-Dimensional Vector
        #============================================#

        #============================================#
        # Gather domain information
        #============================================#

        i_domain_x = (minimum(points_3d[1,:]), maximum(points_3d[1,:]))
        i_domain_y = (minimum(points_3d[2,:]), maximum(points_3d[2,:]))
        i_domain_z = (minimum(points_3d[3,:]), maximum(points_3d[3,:]))

        println("Domain is $i_domain_x × $i_domain_y × $i_domain_z.")

        n_samples_x = ceil(Int, (i_domain_x[end] - i_domain_x[1]) / sampling_rate[1])
        n_samples_y = ceil(Int, (i_domain_y[end] - i_domain_y[1]) / sampling_rate[2])
        n_samples_z = ceil(Int, (i_domain_z[end] - i_domain_z[1]) / sampling_rate[3])

        println("Sampling into $n_samples_x × $n_samples_y × $n_samples_z cells.")

        #============================================#
        # Calculate problem size & hardware limitations
        #============================================#

        n_threads                   = nthreads()

        b_mem_per_var_and_timestep  = sizeof(Float64) * n_samples_x * n_samples_y
        b_mem_static                = sizeof(UInt16)  * n_samples_x * n_samples_y + 1024^2*n_threads

        n_vars_3d                   = length(var_names_3d)
        n_timesteps_per_iteration   = (mem_limit - b_mem_static) ÷ (n_vars_3d * b_mem_per_var_and_timestep)
        n_iterations                = ceil(Int, length(times) / n_timesteps_per_iteration)
        
        n_tetrahedra                = size(tetrahedra, 2)
        n_tetrahedra_per_thread     = n_tetrahedra / n_threads
        n_total_problems            = length(times) * n_tetrahedra

        #============================================#
        # Put tetrahedra into bins to avoid thread
        # synchronization & lock contention later on.
        #============================================#

        println("Binning $n_tetrahedra tetrahedra into $(n_threads)+$(n_threads - 1)+1 buckets...")

        # For each tetrahedron, assign a bin with an ID ∈ 1:2*n_threads based on its location
        l_bin_ids = Array{UInt8, 1}(undef, n_tetrahedra)

        Threads.@threads for thread_id ∈ 1:n_threads
            thread_range_min = floor(INDEX_TYPE, (thread_id - 1) * n_tetrahedra_per_thread) + 1
            thread_range_max = thread_id == n_threads ? n_tetrahedra : floor(INDEX_TYPE, thread_id * n_tetrahedra_per_thread)

            m43_tet_points = Array{Float64, 2}(undef, (4, 3))
            m32_tet_aabb   = Array{Float64, 2}(undef, (3, 2))

            for tet_id ∈ thread_range_min:thread_range_max
                tetrahedron = (@view tetrahedra[:, tet_id])
                for i ∈ 1:4, dim ∈ 1:3
                    m43_tet_points[i, dim] = points_3d[dim, tetrahedron[i]]
                end

                coord_x_min = minimum(m43_tet_points[:, 1])
                coord_x_max = maximum(m43_tet_points[:, 1])

                @inbounds idx_x_min = floor(Int,(coord_x_min - i_domain_x[1]) / sampling_rate[1]) + 1
                @inbounds idx_x_max =  ceil(Int,(coord_x_max - i_domain_x[1]) / sampling_rate[1])

                bin_id_l = floor(UInt8, idx_x_min / n_samples_x * n_threads) + 1
                bin_id_r =  ceil(UInt8, idx_x_max / n_samples_x * n_threads)

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
            end
        end

        # The number of tetrahedra in each respective bin
        l_bin_counts = zeros(INDEX_TYPE, 2*n_threads)
        for bin_id ∈ l_bin_ids
            l_bin_counts[bin_id] += 1
        end

        for bin_id ∈ 1:length(l_bin_counts)
            println(l_bin_counts[bin_id], " Tetrahedra in bin ", bin_id)
        end

        println("Done.")

        #============================================#
        # Set up NetCDF output file
        #============================================#

        nccreate(out_filename, var_names_3d[1], 
                 "y", [i_domain_y[1] + i * sampling_rate[2] for i ∈ 1:n_samples_y], 
                 "x", [i_domain_x[1] + i * sampling_rate[1] for i ∈ 1:n_samples_x], 
                 "time", times, Dict("units"=>"seconds"))

        for var_id ∈ 2:n_vars_3d
            nccreate(out_filename, var_names_3d[var_id], 
                     "y", "x", "time")
        end

        # These are the grids that will be written to the NetCDF output later on.
        l_iter_output_grids = Array{Array{Float64, 3}, 1}(undef, n_vars_3d)
        for var_id ∈ 1:n_vars_3d
            l_iter_output_grids[var_id] = Array{Float64, 3}(undef, (n_samples_y, n_samples_x, n_timesteps_per_iteration))
        end

        # This matrix counts the total number of samples in each output grid cell
        # (read: the total number of sampled points in z-direction for each xy-index of the output grids)
        # This stays constant across all vars and timesteps because the grid geometry is constant.
        myx_output_sample_counts = Array{UInt16, 2}(undef, (n_samples_y, n_samples_x))

        #============================================#
        # Sample tetrahedra and write to output
        #============================================#

        println("RAM limit: ", Main.Util.human_readable_size(mem_limit), "; Can work on ", n_timesteps_per_iteration, 
                " timesteps at once (therefore needing ", n_iterations , " iterations).")

        # In each iteration, process ALL tetrahedra for ALL variables for as many timesteps as possible*
        # *possible means that the total memory footprint cannot be exceeded. This limits the number of 
        # timesteps that can be stored in l_iter_output_grids. 
        for iteration ∈ 1:n_iterations
            t_start = (iteration - 1) * n_timesteps_per_iteration + 1
            n_times = min(length(times) - (iteration - 1) * n_timesteps_per_iteration, n_timesteps_per_iteration)

            #============================================#
            # Prefetch the values of all vars in all
            # timesteps for the current iteration.
            # This eliminates random accesses later on
            # and D R A S T I C A L L Y improves runtime.
            #============================================#

            prefetched_vars = Array{Float64, 3}(undef, (n_tetrahedra, n_times, n_vars_3d))
            @profile begin
                for var_id ∈ 1:n_vars_3d, time ∈ t_start:t_start + n_times - 1
                    dest_offset = (var_id-1)*n_times*n_tetrahedra + (time-t_start)*n_tetrahedra + 1
                    copyto!(prefetched_vars, dest_offset, vars_3d[time, var_id], 1, n_tetrahedra)
                end
            end # profile

            # Reset the output values to 0. Only do so for the timesteps actually processed in this iteration.
            for grid in l_iter_output_grids
                grid[1:n_times,:,:] .= 0
            end

            #============================================#
            # Binning
            #============================================#
            # First, process each of the n_threads normal
            # bins.
            #
            # Then, process the n_threads-1 bins on the
            # borders between normal bins.
            #
            # Lastly, process the 1 bin of tetrahedra
            # that overlapped more than 2 adjacent bins.
            #============================================#

            Threads.@threads for thread_id ∈ 1:n_threads
                iterate_tets!(
                    thread_id,         l_bin_ids,   l_bin_counts,
                    n_tetrahedra,      tetrahedra,  points_3d, prefetched_vars,
                    i_domain_x,        i_domain_y,  i_domain_z,
                    n_samples_x,       n_samples_y, n_samples_z,
                    sampling_rate,     t_start,     n_times,
                    l_iter_output_grids, myx_output_sample_counts, print_progress = (thread_id == n_threads ÷ 2))
            end # for thread_id

            println("Processing tets on bucket borders...")

            # Now, iterate over the remaining tets that lie ON bin edges and therefore could not be processed in parallel
            # with the normal bins without conflict
            Threads.@threads for thread_id ∈ 1:n_threads-1
                iterate_tets!(
                    n_threads + thread_id, l_bin_ids,   l_bin_counts,
                    n_tetrahedra,          tetrahedra,  points_3d, prefetched_vars,
                    i_domain_x,            i_domain_y,  i_domain_z,
                    n_samples_x,           n_samples_y, n_samples_z,
                    sampling_rate,         t_start,     n_times,
                    l_iter_output_grids,     myx_output_sample_counts, print_progress = (thread_id == n_threads ÷ 2))
            end # for thread_id

            # Lastly, the tets that overlapped multiple bins (practically never occurs!)
            iterate_tets!(
                2 * n_threads,     l_bin_ids,   l_bin_counts,
                n_tetrahedra,      tetrahedra,  points_3d, prefetched_vars,
                i_domain_x,        i_domain_y,  i_domain_z,
                n_samples_x,       n_samples_y, n_samples_z,
                sampling_rate,     t_start,     n_times,
                l_iter_output_grids, myx_output_sample_counts, print_progress = true)

            println("Done.")

            #============================================#
            # Write NetCDF outputs for this iteration
            #============================================#

            println("Writing outputs for timesteps $(t_start-1)-$(t_start+n_times-2)...")
            
            start = [1, 1, t_start]
            count = [-1, -1, n_times]

            for var_id ∈ 1:n_vars_3d
                ncwrite(l_iter_output_grids[var_id], out_filename, var_names_3d[var_id], start=start, count=count)
            end

            println("Done.")
        end # for iteration

    end # function rasterize

    function iterate_tets!(
        bin_id,             l_bin_ids,      l_bin_counts, 
        n_tetrahedra,       l_tetrahedra,   l_points, vars, 
        i_domain_x,         i_domain_y,     i_domain_z,
        n_samples_x,        n_samples_y,    n_samples_z,
        v_sampling_rate,    t_start,        n_times,
        l_iteration_grids,  myx_grid_sample_counts; print_progress=false)

        m43_tet_points      = Array{Float64, 2}(undef, (4, 3))
        m32_tet_aabb        = Array{Float64, 2}(undef, (3, 2))
        tet_sample_counts   = Array{UInt8, 2}(undef, (64, 64))
        
        l_face_points       = Array{SubArray, 1}(undef, 4)
        l_tet_faces         = Array{TetrahedronFace, 1}(undef, 4)

        v_p                 = Array{Float64, 1}(undef, 3)
        v_n                 = Array{Float64, 1}(undef, 3)

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

            l_tetrahedron = (@view l_tetrahedra[:, tet_id])
            for i ∈ 1:4, dim ∈ 1:3
                m43_tet_points[i, dim] = l_points[dim, l_tetrahedron[i]]
            end

            for dim ∈ 1:3
                m32_tet_aabb[dim, 1] = minimum(m43_tet_points[:, dim])
                m32_tet_aabb[dim, 2] = maximum(m43_tet_points[:, dim])
            end

            #============================================#
            # A note on global & local indices
            #============================================#
            # In the following, idx_g_ will refer to an
            # index of the global sampling domain while
            # idx_l_ refers to an index in the tet-local
            # domain.
            #============================================#

            @inbounds idx_g_x_min = floor(Int32, (m32_tet_aabb[1,1] - i_domain_x[1]) / v_sampling_rate[1]) + 1
            @inbounds idx_g_y_min = floor(Int32, (m32_tet_aabb[2,1] - i_domain_y[1]) / v_sampling_rate[2]) + 1
            @inbounds idx_g_z_min = floor(Int32, (m32_tet_aabb[3,1] - i_domain_z[1]) / v_sampling_rate[3]) + 1
     
            @inbounds idx_g_x_max =  ceil(Int32, (m32_tet_aabb[1,2] - i_domain_x[1]) / v_sampling_rate[1])
            @inbounds idx_g_y_max =  ceil(Int32, (m32_tet_aabb[2,2] - i_domain_y[1]) / v_sampling_rate[2])
            @inbounds idx_g_z_max =  ceil(Int32, (m32_tet_aabb[3,2] - i_domain_z[1]) / v_sampling_rate[3])

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

                    @printf("Working on tetrahedron %d of %d in timesteps %d-%d (%2.2f%% done). ETR: %02d:%02d:%02d\n", 
                            n_bin_rasterized, n_bin_total,
                            t_start - 1,           # Tools like ParaView work 0-indexed. Avoid confusion by outputting 0-indexed here.
                            t_start + n_times - 2, # Tools like ParaView work 0-indexed. Avoid confusion by outputting 0-indexed here.
                            n_bin_rasterized/n_bin_total*100,
                            hh.value, mm.value, ss.value)
                end
            end # if thread_id

            n_current_cells_x = idx_g_x_max - idx_g_x_min + 1
            n_current_cells_y = idx_g_y_max - idx_g_y_min + 1

            # If the current tet has a bounding box larger than anticipated, enlarge the sample counts array accordingly
            # By allocating double the size of the bounding box, future tets that are just a bit larger than the current one
            # will not trigger a resize operation.
            if size(tet_sample_counts, 1) < n_current_cells_y || size(tet_sample_counts, 2) < n_current_cells_x
                tet_sample_counts = Array{UInt8, 2}(undef, (n_current_cells_y * 2, n_current_cells_x * 2))
            end

            # For the bounding box of the current tet, reset all sample counts to 0
            @inbounds tet_sample_counts[1:n_current_cells_y, 1:n_current_cells_x] .= 0

            for i ∈ 1:4
                @inbounds l_face_points[i] = @view m43_tet_points[i,:]
            end

            # Array{(n, d, dist_p_excl)}
            for i ∈ 1:4
                # Exclude point i from the tetrahedron and consider the plane defined by the other 3 points
                # If and only if (x,y,z) lies on the "positive" side (in normal direction) of all 4 planes, 
                # it is inside the tetrahedron
                @inbounds v_p_excl =    l_face_points[i]
                @inbounds v_p1 =        l_face_points[(i%4)+1]
                @inbounds v_p2 =        l_face_points[((i+1)%4)+1]
                @inbounds v_p3 =        l_face_points[((i+2)%4)+1]
            
                # Calculate Hesse normal form of the plane defined by p1, p2, p3
                # Normal vector
                cross3!(v_p2-v_p1, v_p3-v_p1, v_n)
                normalize!(v_n)
            
                # n1x1 + n2x2 + n3x3 + d = 0; solve for d
                d = -dot3(v_p1, v_n)
                dist_p_excl = dot3(v_p_excl, v_n) + d

                # Ensure that p_excl is considered to be on the "positive" side of the plane
                if dist_p_excl < 0
                    v_n = -v_n
                    d = -d
                    dist_p_excl = -dist_p_excl
                end

                @inbounds l_tet_faces[i] = TetrahedronFace((v_n[1], v_n[2], v_n[3]), d, dist_p_excl)
            end # for i

            #============================================#
            # Sample each cuboid that is overlapping the 
            # current tet's AABB.
            # If its CENTER POINT lies within the tet,
            # increase the tet's sample count for the
            # cuboid's (x, y) index by 1.
            #============================================#

            for idx_g_x ∈ idx_g_x_min:idx_g_x_max
                @inbounds v_p[1] = (idx_g_x-1)*v_sampling_rate[1] + i_domain_x[1]
                for idx_g_y ∈ idx_g_y_min:idx_g_y_max
                    @inbounds v_p[2] = (idx_g_y-1)*v_sampling_rate[2] + i_domain_y[1]
                    for idx_g_z ∈ idx_g_z_min:idx_g_z_max
                        @inbounds v_p[3] = (idx_g_z-1)*v_sampling_rate[3] + i_domain_z[1]
                    
                        accepted = true
                        # For each tetrahedron face, filter out cells that are not in the tetrahedron
                        for i ∈ 1:4
                            @inbounds face = l_tet_faces[i]
                            # p_excl is ALWAYS on the "positive" side of the plane. 
                            # So, if p_excl and p lie on the same side, p is inside the tetrahedron
                            dist_p = dot3(v_p, face.n) + face.d
                        
                            # The point lies approx. ON the face, reject in the tet on the lexicographically
                            # more negative side of the face. This ensures that for two neighboring tetrahedra, the point is
                            # only accepted for one of them.
                            # Def.: lexicographical order: 
                            #       a ≤ b ⟺ (a.x ≤ b.x) ∨ (a.x = b.x ∧ a.y ≤ b.y) ∨ (a.x = b.x ∧ a.y = b.y ∧ a.z ≤ b.z)
                            # This is a total ordering on all vectors in R^3 and will therefore always be larger for the "positive"
                            # normal vector in contrast to its negated counterpart.
                            if (abs(dist_p) <= EDGE_TOL) 
                                @inbounds is_positive_side = 
                                    (v_n[1] > 0.) || 
                                    (v_n[1] == 0. && v_n[2] > 0.) || 
                                    (v_n[1] == 0. && v_n[2] == 0. && v_n[3] > 0.)

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

                        if accepted # accept cell and increment its 2D cell counter for the current tetrahedron by 1
                            idx_l_x = idx_g_x - idx_g_x_min + 1
                            idx_l_y = idx_g_y - idx_g_y_min + 1
                            @inbounds tet_sample_counts[idx_l_y, idx_l_x] += 1
                        end
                    end # for z
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

            for idx_l_x ∈ 1:n_current_cells_x, idx_l_y ∈ 1:n_current_cells_y
                @inbounds num_samples = tet_sample_counts[idx_l_y, idx_l_x]
                if num_samples > 0
                    idx_g_x = idx_l_x + idx_g_x_min - 1; idx_g_y = idx_l_y + idx_g_y_min - 1
                    @inbounds total_samples = (myx_grid_sample_counts[idx_g_y, idx_g_x] += num_samples)
                    for var_id ∈ 1:3
                        for t ∈ 1:n_times
                            @inbounds avg = l_iteration_grids[var_id][idx_g_y, idx_g_x, t]
                            @inbounds l_iteration_grids[var_id][idx_g_y, idx_g_x, t] = avg + (vars[tet_id, t, var_id]-avg)*(num_samples/total_samples)
                        end
                    end
                end
            end

            n_bin_rasterized += 1
        end # for tet_id
    end

    @inline function cross3!(a, b, ret)
        @inbounds ret[1] = (a[2]*b[3] - a[3]*b[2])
        @inbounds ret[2] = (a[3]*b[1] - a[1]*b[3])
        @inbounds ret[3] = (a[1]*b[2] - a[2]*b[1])
    end

    @inline function dot3(a, b)
        @inbounds return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
    end
end # module Rasterization
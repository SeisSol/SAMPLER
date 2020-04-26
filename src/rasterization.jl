module Rasterization
    using Base.Threads
    using LinearAlgebra
    using Dates
    using Printf
    using NetCDF
    using Profile
    using Traceur
    using Main.Util

    struct HesseNormalForm
        v_n0 :: NTuple{3, Float64}
        d    :: Float64
    end

    mutable struct ThreadContext{DIMS}
        mn3_simplex_points    :: Array{Float64, 2}
        m32_simplex_aabb      :: Array{Float64, 2}
        simplex_sample_counts :: Array{UInt8, 2}
        l_face_points         :: Array{SubArray, 1}
        l_simplex_faces       :: Array{HesseNormalForm, 1}
        v_p                   :: Array{Float64, 1}
        v_n                   :: Array{Float64, 1}

        d_start               :: DateTime
        d_last_printed        :: DateTime
        print_interval        :: Second
        n_bin_rasterized      :: Int
        n_bin_total           :: Int

        thread_range_min      :: Int
        thread_range_max      :: Int

        function ThreadContext{DIMS}() where DIMS
            mn3_simplex_points    = Array{Float64, 2}(undef, (DIMS + 1, 3))
            m32_simplex_aabb      = Array{Float64, 2}(undef, (3, 2))
            simplex_sample_counts = Array{UInt8, 2}(undef, (64, 64))
            l_face_points         = Array{SubArray, 1}(undef, DIMS + 1)
            l_simplex_faces       = Array{HesseNormalForm, 1}(undef, DIMS + 1)
            v_p                   = Array{Float64, 1}(undef, 3)
            v_n                   = Array{Float64, 1}(undef, 3)

            new(
                mn3_simplex_points    ,
                m32_simplex_aabb      ,
                simplex_sample_counts ,
                l_face_points         ,
                l_simplex_faces       ,
                v_p                   ,
                v_n
            ) 
        end
    end

    mutable struct RasterizationContext{DIMS, INDEX_TYPE <: Integer}
        dims                       :: Int

        vn_i_domain                :: NTuple{DIMS, NTuple{2, Float64}}
        vn_n_samples               :: NTuple{DIMS,Int}
        v_sampling_rate            :: NTuple{3, Float64}

        n_threads                  :: Int

        b_mem_per_var_and_timestep :: Int
        b_mem_static               :: Int

        n_vars                     :: Int
        mtv_vars                   :: AbstractArray{A where A <: AbstractArray, 2}

        n_timesteps_per_iteration  :: Int
        n_iterations               :: Int
        
        n_simplices                :: Int
        mnt_simplices              :: AbstractArray{INDEX_TYPE, 2}
        m3t_points                 :: AbstractArray{Float64, 2}

        simplices_per_thread       :: Float64
        n_total_problems           :: Int

        l_bin_ids                  :: Array{UInt8, 1}
        l_bin_counts               :: Array{Int32, 1}

        l_iter_output_grids        :: Array{Array{Float64, 3}, 1}
        myx_output_sample_counts   :: Array{UInt16, 2}

        l_thread_contexts          :: Array{ThreadContext{DIMS}, 1}

        function RasterizationContext{DIMS, INDEX_TYPE}(
            simplices, points, vars, var_names, times, sampling_rate, mem_limit) where DIMS where INDEX_TYPE <: Integer

            mtv_vars        = vars
            mnt_simplices   = simplices
            m3t_points      = points
        
            #============================================#
            # Gather domain information
            #============================================#

            v_sampling_rate = sampling_rate
            vn_i_domain     = Tuple((minimum(points[dim,:]), maximum(points[dim,:])) for dim ∈ 1:DIMS)
            vn_n_samples    = Tuple(ceil(Int, (vn_i_domain[dim][end] - vn_i_domain[dim][1]) / v_sampling_rate[dim]) for dim ∈ 1:3)

            #============================================#
            # Calculate problem size & hardware limitations
            #============================================#

            n_threads       = nthreads()

            b_mem_per_var_and_timestep = sizeof(Float64) * vn_n_samples[1] * vn_n_samples[2]
            b_mem_static               = sizeof(UInt16)  * vn_n_samples[1] * vn_n_samples[2] + 1024^2 * n_threads
        
            n_vars                     = length(var_names)
            n_timesteps_per_iteration  = (mem_limit - b_mem_static) ÷ (n_vars * b_mem_per_var_and_timestep)
            n_iterations               = ceil(Int, length(times) / n_timesteps_per_iteration)
        
            n_simplices                = size(simplices, 2)
            simplices_per_thread       = n_simplices / n_threads
            n_total_problems           = length(times) * n_simplices

            #============================================#
            # Allocate big, shared arrays for later use
            #============================================#

            l_bin_ids    = Array{UInt8, 1}(undef, n_simplices)
            l_bin_counts = zeros(INDEX_TYPE, 2 * n_threads)

            l_thread_contexts = Array{ThreadContext{DIMS}, 1}(undef, n_threads)

            l_iter_output_grids      = Array{Array{Float64, 3}, 1}(undef, n_vars)
            myx_output_sample_counts = Array{UInt16, 2}(undef, (vn_n_samples[2], vn_n_samples[1]))

            new(DIMS                       ,
                vn_i_domain                ,
                vn_n_samples               ,
                v_sampling_rate            ,
                n_threads                  ,
                b_mem_per_var_and_timestep ,
                b_mem_static               ,
                n_vars                     ,
                mtv_vars                   ,
                n_timesteps_per_iteration  ,
                n_iterations               ,
                n_simplices                ,
                mnt_simplices              ,
                m3t_points                 ,
                simplices_per_thread       ,
                n_total_problems           ,
                l_bin_ids                  ,
                l_bin_counts               ,
                l_iter_output_grids        ,
                myx_output_sample_counts   ,
                l_thread_contexts          
            )
        end
    end

    mutable struct IterationContext
        t_start         :: Int
        n_times         :: Int
        prefetched_vars :: Array{Float64, 3}
        nc_start        :: Array{Int, 1}
        nc_count        :: Array{Int, 1}

        function IterationContext(cr :: RasterizationContext)
            prefetched_vars = Array{Float64, 3}(undef, (cr.n_simplices, cr.n_timesteps_per_iteration, cr.n_vars))
            nc_start = [1, 1, 0]
            nc_count = [-1, -1, 0]

            new(0, 0, prefetched_vars, nc_start, nc_count)
        end
    end

    const EDGE_TOL = .0

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
        #   t_  : Discrete Timestep (Integer)
        #   cr  : Rasterization Context
        #   ct  : Thread Context
        #   ci  : Iteration Context
        #============================================#

        cr = RasterizationContext{3, INDEX_TYPE}(tetrahedra, points_3d, vars_3d, var_names_3d, times, sampling_rate, mem_limit)

        println("Domain is $(cr.vn_i_domain[1]) × $(cr.vn_i_domain[2]) × $(cr.vn_i_domain[3]).")
        println("Sampling into $(cr.vn_n_samples[1]) × $(cr.vn_n_samples[2]) × $(cr.vn_n_samples[3]) cells.")

        #============================================#
        # Put tetrahedra into bins to avoid thread
        # synchronization & lock contention later on.
        #============================================#

        println("Binning $(cr.n_simplices) tetrahedra into $(cr.n_threads)+$(cr.n_threads - 1)+1 buckets...")

        # For each tetrahedron, assign a bin with an ID ∈ 1:2*n_threads based on its location
        Threads.@threads for thread_id ∈ 1:cr.n_threads
            ct = cr.l_thread_contexts[thread_id] = ThreadContext{cr.dims}()

            ct.thread_range_min = floor(INDEX_TYPE, (thread_id - 1) * cr.simplices_per_thread) + 1
            ct.thread_range_max = thread_id == cr.n_threads ? cr.n_simplices : floor(INDEX_TYPE, thread_id * cr.simplices_per_thread)

            for tet_id ∈ ct.thread_range_min:ct.thread_range_max
                tetrahedron = (@view cr.mnt_simplices[:, tet_id])
                for i ∈ 1:(cr.dims + 1), dim ∈ 1:3
                    ct.mn3_simplex_points[i, dim] = cr.m3t_points[dim, tetrahedron[i]]
                end

                cr.l_bin_ids[tet_id] = bin_of(ct.mn3_simplex_points, cr)
            end
        end

        # The number of tetrahedra in each respective bin
        for bin_id ∈ cr.l_bin_ids
            cr.l_bin_counts[bin_id] += 1
        end

        for bin_id ∈ 1:length(cr.l_bin_counts)
            println(cr.l_bin_counts[bin_id], " Tetrahedra in bin ", bin_id)
        end

        println("Done.")

        #============================================#
        # Set up NetCDF output file
        #============================================#

        nccreate(out_filename, var_names_3d[1], 
                 "y", [cr.vn_i_domain[2][1] + i * cr.v_sampling_rate[2] for i ∈ 1:cr.vn_n_samples[2]], 
                 "x", [cr.vn_i_domain[1][1] + i * cr.v_sampling_rate[1] for i ∈ 1:cr.vn_n_samples[1]], 
                 "time", times, Dict("units"=>"seconds"))

        for var_name ∈ var_names_3d[2:end]
            nccreate(out_filename, var_name, "y", "x", "time")
        end

        # These are the grids that will be written to the NetCDF output later on.
        
        for var_id ∈ 1:cr.n_vars
            cr.l_iter_output_grids[var_id] = Array{Float64, 3}(undef, (cr.vn_n_samples[2], cr.vn_n_samples[1], 
                                                               cr.n_timesteps_per_iteration))
        end

        # This matrix counts the total number of samples in each output grid cell
        # (read: the total number of sampled points in z-direction for each xy-index of the output grids)
        # This stays constant across all vars and timesteps because the grid geometry is constant.

        #============================================#
        # Sample tetrahedra and write to output
        #============================================#

        println("RAM limit: ", Main.Util.human_readable_size(mem_limit), "; Can work on ", cr.n_timesteps_per_iteration, 
                " timesteps at once (therefore needing ", cr.n_iterations , " iterations).")

        # In each iteration, process ALL tetrahedra for ALL variables for as many timesteps as possible*
        # *possible means that the total memory footprint cannot be exceeded. This limits the number of 
        # timesteps that can be stored in l_iter_output_grids. 
        ci = IterationContext(cr)

        for iteration ∈ 1:cr.n_iterations

            ci.t_start = (iteration - 1) * cr.n_timesteps_per_iteration + 1
            ci.n_times = min(length(times) - (iteration - 1) * cr.n_timesteps_per_iteration, cr.n_timesteps_per_iteration)

            #============================================#
            # Prefetch the values of all vars in all
            # timesteps for the current iteration.
            # This eliminates random accesses later on
            # and D R A S T I C A L L Y improves runtime.
            #============================================#

            begin
                for var_id ∈ 1:cr.n_vars, time ∈ ci.t_start:(ci.t_start + ci.n_times - 1)
                    dest_offset = (var_id-1)*ci.n_times*cr.n_simplices + (time-ci.t_start)*cr.n_simplices + 1
                    copyto!(ci.prefetched_vars, dest_offset, cr.mtv_vars[time, var_id], 1, cr.n_simplices)
                end
            end # profile

            # Reset the output values to 0. Only do so for the timesteps actually processed in this iteration.
            for grid in cr.l_iter_output_grids
                grid[1:ci.n_times,:,:] .= 0
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

            Threads.@threads for thread_id ∈ 1:cr.n_threads
                iterate_tets!(thread_id, cr, ci, cr.l_thread_contexts[thread_id], 
                              print_progress = (thread_id == 1))
            end # for thread_id

            open("profile.txt", "w") do p
                Profile.print(p)
            end

            println("Processing tets on bucket borders...")

            # Now, iterate over the remaining tets that lie ON bin edges and therefore could not be processed in parallel
            # with the normal bins without conflict
            Threads.@threads for thread_id ∈ 1:cr.n_threads-1
                iterate_tets!(cr.n_threads + thread_id, cr, ci, cr.l_thread_contexts[thread_id], 
                              print_progress = (thread_id == 1))
            end # for thread_id

            # Lastly, the tets that overlapped multiple bins (practically never occurs!)
            iterate_tets!(2 * cr.n_threads, cr, ci, cr.l_thread_contexts[1], print_progress = true)

            println("Done.")

            #============================================#
            # Write NetCDF outputs for this iteration
            #============================================#

            println("Writing outputs for timesteps $(ci.t_start-1)-$(ci.t_start+ci.n_times-2)...")
            
            ci.nc_start[3] = ci.t_start
            ci.nc_count[3] = ci.n_times

            for var_id ∈ 1:cr.n_vars
                ncwrite(cr.l_iter_output_grids[var_id], out_filename, var_names_3d[var_id], start=ci.nc_start, count=ci.nc_count)
            end

            println("Done.")
        end # for iteration

    end # function rasterize

    function iterate_tets!(bin_id, cr :: RasterizationContext{DIMS, I}, ci :: IterationContext, ct :: ThreadContext{DIMS}; print_progress=false) where DIMS where I
        
        ct.d_start          = now()
        ct.d_last_printed   = ct.d_start
        ct.print_interval   = Second(2)
        ct.n_bin_rasterized = 0
        ct.n_bin_total      = cr.l_bin_counts[bin_id]

        for tet_id ∈ 1:cr.n_simplices
            # Only process tetrahedra in thread's own bin
            if cr.l_bin_ids[tet_id] != bin_id
                continue
            end

            if print_progress
                d_now = now()

                if (d_now - ct.d_last_printed) >= ct.print_interval && ct.n_bin_rasterized > 0
                    ct.d_last_printed = d_now

                    etr  = (d_now-ct.d_start) * (ct.n_bin_total - ct.n_bin_rasterized) ÷ (ct.n_bin_rasterized)
                    hh   = floor(etr, Dates.Hour)
                    etr -= floor(hh, Dates.Millisecond)
                    mm   = floor(etr, Dates.Minute)
                    etr -= floor(mm, Dates.Millisecond)
                    ss   = floor(etr, Dates.Second)

                    @printf("Working on tetrahedron %d of %d in timesteps %d-%d (%2.2f%% done). ETR: %02d:%02d:%02d\n", 
                            ct.n_bin_rasterized, ct.n_bin_total,
                            ci.t_start - 1,              # Tools like ParaView work 0-indexed
                            ci.t_start + ci.n_times - 2, # ⇒ Avoid confusion by outputting 0-indexed here.
                            ct.n_bin_rasterized/ct.n_bin_total*100,
                            hh.value, mm.value, ss.value)
                end
            end # if print_progress

            n_current_cells_x, n_current_cells_y, idx_g_x_min, idx_g_y_min = rasterize_simplex(cr, ci, ct, tet_id)

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
                #=@inbounds=# num_samples = ct.simplex_sample_counts[idx_l_y, idx_l_x]
                if num_samples > 0
                    idx_g_x = idx_l_x + idx_g_x_min - 1; idx_g_y = idx_l_y + idx_g_y_min - 1
                    #=@inbounds=# total_samples = (cr.myx_output_sample_counts[idx_g_y, idx_g_x] += num_samples)
                    for var_id ∈ 1:3
                        for t ∈ 1:ci.n_times
                            #=@inbounds=# avg = cr.l_iter_output_grids[var_id][idx_g_y, idx_g_x, t]
                            #=@inbounds=# cr.l_iter_output_grids[var_id][idx_g_y, idx_g_x, t] = avg + (ci.prefetched_vars[tet_id, t, var_id]-avg)*(num_samples/total_samples)
                        end
                    end
                end
            end

            ct.n_bin_rasterized += 1
        end # for tet_id
        println("Thread $(threadid()) done.")
    end

    function rasterize_simplex(cr :: RasterizationContext, ci :: IterationContext, ct :: ThreadContext{DIMS}, simplex_id :: Int) where DIMS
        #=@inbounds=# l_simplex = (@view cr.mnt_simplices[:, simplex_id])
        n_points = DIMS + 1

        for i ∈ 1:n_points, dim ∈ 1:3
            #=@inbounds=# ct.mn3_simplex_points[i, dim] = cr.m3t_points[dim, l_simplex[i]]
        end

        for dim ∈ 1:3
            #=@inbounds=# ct.m32_simplex_aabb[dim, 1] = minimum(ct.mn3_simplex_points[:, dim])
            #=@inbounds=# ct.m32_simplex_aabb[dim, 2] = maximum(ct.mn3_simplex_points[:, dim])
        end

        #============================================#
        # A note on global & local indices
        #============================================#
        # In the following, idx_g_ will refer to an
        # index of the global sampling domain while
        # idx_l_ refers to an index in the tet-local
        # domain.
        #============================================#

        #=@inbounds=# idx_g_x_min = floor(Int32, (ct.m32_simplex_aabb[1,1] - cr.vn_i_domain[1][1]) / cr.v_sampling_rate[1]) + 1
        #=@inbounds=# idx_g_y_min = floor(Int32, (ct.m32_simplex_aabb[2,1] - cr.vn_i_domain[2][1]) / cr.v_sampling_rate[2]) + 1
        #=@inbounds=# idx_g_z_min = floor(Int32, (ct.m32_simplex_aabb[3,1] - cr.vn_i_domain[3][1]) / cr.v_sampling_rate[3]) + 1
    
        #=@inbounds=# idx_g_x_max =  ceil(Int32, (ct.m32_simplex_aabb[1,2] - cr.vn_i_domain[1][1]) / cr.v_sampling_rate[1])
        #=@inbounds=# idx_g_y_max =  ceil(Int32, (ct.m32_simplex_aabb[2,2] - cr.vn_i_domain[2][1]) / cr.v_sampling_rate[2])
        #=@inbounds=# idx_g_z_max =  ceil(Int32, (ct.m32_simplex_aabb[3,2] - cr.vn_i_domain[3][1]) / cr.v_sampling_rate[3])

        n_current_cells_x = idx_g_x_max - idx_g_x_min + 1
        n_current_cells_y = idx_g_y_max - idx_g_y_min + 1

        # If the current tet has a bounding box larger than anticipated, enlarge the sample counts array accordingly
        # By allocating double the size of the bounding box, future tets that are just a bit larger than the current one
        # will not trigger a resize operation.
        if size(ct.simplex_sample_counts, 1) < n_current_cells_y || size(ct.simplex_sample_counts, 2) < n_current_cells_x
            ct.simplex_sample_counts = Array{UInt8, 2}(undef, (n_current_cells_y * 2, n_current_cells_x * 2))
        end

        # For the bounding box of the current tet, reset all sample counts to 0
        #=@inbounds=# ct.simplex_sample_counts[1:n_current_cells_y, 1:n_current_cells_x] .= 0

        for i ∈ 1:n_points
            #=@inbounds=# ct.l_face_points[i] = @view ct.mn3_simplex_points[i,:]
        end

        calculate_hnfs!(ct)

        #============================================#
        # Sample each cuboid that is overlapping the 
        # current tet's AABB.
        # If its CENTER POINT lies within the tet,
        # increase the tet's sample count for the
        # cuboid's (x, y) index by 1.
        #============================================#

        for idx_g_x ∈ idx_g_x_min:idx_g_x_max
            #=@inbounds=# ct.v_p[1] = (idx_g_x-1)*cr.v_sampling_rate[1] + cr.vn_i_domain[1][1]

            for idx_g_y ∈ idx_g_y_min:idx_g_y_max
                #=@inbounds=# ct.v_p[2] = (idx_g_y-1)*cr.v_sampling_rate[2] + cr.vn_i_domain[2][1]

                for idx_g_z ∈ idx_g_z_min:idx_g_z_max
                    #=@inbounds=# ct.v_p[3] = (idx_g_z-1)*cr.v_sampling_rate[3] + cr.vn_i_domain[3][1]
                
                    accepted = true
                    # For each tetrahedron face, filter out cells that are not in the tetrahedron
                    for i ∈ 1:n_points
                        #=@inbounds=# face = ct.l_simplex_faces[i]
                        # p_excl is ALWAYS on the "positive" side of the plane. 
                        # So, if p_excl and p lie on the same side, p is inside the tetrahedron
                        dist_p = dot3(ct.v_p, face.v_n0) + face.d
                    
                        # The point lies approx. ON the face, reject in the tet on the lexicographically
                        # more negative side of the face. This ensures that for two neighboring tetrahedra, the point is
                        # only accepted for one of them.
                        # Def.: lexicographical order: 
                        #       a ≤ b ⟺ (a.x ≤ b.x) ∨ (a.x = b.x ∧ a.y ≤ b.y) ∨ (a.x = b.x ∧ a.y = b.y ∧ a.z ≤ b.z)
                        # This is a total ordering on all vectors in R^3 and will therefore always be larger for the "positive"
                        # normal vector in contrast to its negated counterpart.
                        if (abs(dist_p) <= EDGE_TOL) 
                            #=@inbounds=# is_positive_side = 
                                (face.v_n0[1] > 0.) || 
                                (face.v_n0[1] == 0. && face.v_n0[2] > 0.) || 
                                (face.v_n0[1] == 0. && face.v_n0[2] == 0. && face.v_n0[3] > 0.)

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
                        #=@inbounds=# ct.simplex_sample_counts[idx_l_y, idx_l_x] += 1
                    end
                end # for z
            end # for y
        end # for x

        return n_current_cells_x, n_current_cells_y, idx_g_x_min, idx_g_y_min
    end

    function bin_of(mn3_points :: AbstractArray{Float64, 2}, cr :: RasterizationContext)
        #=@inbounds=# coord_x_min = minimum(mn3_points[:, 1])
        #=@inbounds=# coord_x_max = maximum(mn3_points[:, 1])

        #=@inbounds=# idx_x_min = floor(Int,(coord_x_min - cr.vn_i_domain[1][1]) / cr.v_sampling_rate[1]) + 1
        #=@inbounds=# idx_x_max =  ceil(Int,(coord_x_max - cr.vn_i_domain[1][1]) / cr.v_sampling_rate[1])

        #=@inbounds=# bin_id_l = floor(UInt8, idx_x_min / cr.vn_n_samples[1] * cr.n_threads) + 1
        #=@inbounds=# bin_id_r =  ceil(UInt8, idx_x_max / cr.vn_n_samples[1] * cr.n_threads)

        # n_bins := 2*n_columns
        # bin i ∈ 1:n_columns contains the tets completely inside the column of bin i
        # bin j ∈ n_columns+1:2*n_columns-1 contains tets that are overlapping across bin i and i+1
        # bin k = 2*n_columns contains tets overlapping over more than 2 adjascent columns
        if bin_id_l == bin_id_r
            return bin_id_l
        elseif bin_id_l == bin_id_r - 1
            return cr.n_threads + bin_id_l
        else
            return 2 * cr.n_threads
        end
    end

    function calculate_hnfs!(ct :: ThreadContext{DIMS}) where DIMS
        n_points = DIMS + 1

        for i ∈ 1:n_points
            # Exclude point i from the tetrahedron and consider the plane defined by the other 3 points
            # If and only if (x,y,z) lies on the "positive" side (in normal direction) of all 4 planes, 
            # it is inside the tetrahedron
            #=@inbounds=# v_p_excl = ct.l_face_points[  i                   ]
            #=@inbounds=# v_p1 =     ct.l_face_points[( i    % n_points) + 1]
            #=@inbounds=# v_p2 =     ct.l_face_points[((i+1) % n_points) + 1]
            #=@inbounds=# v_p3 =     
                if DIMS == 3
                    ct.l_face_points[((i+2) % n_points) + 1]
                elseif DIMS == 2
                    # We are in the 2D case. How do we get the normal vector of a line?
                    # The same way we get it from a plane but the 3rd point is directly above the 1st point,
                    # making p1p2 and p1p3 orthogonal to each other and the result of the cross product later
                    # will be a vector orthogonal to p1p2 and p1p3 - done.
                    [v_p1[1], v_p1[2], v_p1[3] + 1.]
                else
                    throw(DimensionMismatch("Dimension $(n_points - 1) invalid. Only 2 (triangles) and 3 (tetrahedra) are allowed."))
                end
        
            # Calculate Hesse Normal Form (HNF) of the plane defined by p1, p2, p3
            #            n0 ∘ (p - p1) = 0
            # ⟺      n0 ∘ p - n0 ∘ p1 = 0
            # ⟺               n0 ∘ p  = n0 ∘ p1 =: d
            # ⟺               n0 ∘ p  = d

            # Normal vector
            cross3!(v_p2-v_p1, v_p3-v_p1, ct.v_n)

            # In the 2D case, the 3rd component of n is set to 0. This should be the case anyways due to the normal vector hack
            # explained above. Just making sure.
            DIMS == 2 && @assert (ct.v_n[3] == 0.) "2D normal vector has a non-zero z component!"

            normalize!(ct.v_n)
        
            # Remember that n3 = 0 if we are in 2D (n_points = 3)
            # n1x1 + n2x2 + (n3x3) + d = 0; solve for d
            d = -dot3(v_p1, ct.v_n)
            dist_p_excl = dot3(v_p_excl, ct.v_n) + d

            # Ensure that p_excl is considered to be on the "positive" side of the plane
            # This is equivalent to the normal vector pointing INTO the simplex
            if dist_p_excl < 0
                ct.v_n = -ct.v_n
                d = -d
            end

            #=@inbounds=# ct.l_simplex_faces[i] = HesseNormalForm((ct.v_n[1], ct.v_n[2], ct.v_n[3]), d)
        end # for i
    end

    @inline function cross3!(a, b, ret)
        #=@inbounds=# ret[1] = (a[2]*b[3] - a[3]*b[2])
        #=@inbounds=# ret[2] = (a[3]*b[1] - a[1]*b[3])
        #=@inbounds=# ret[3] = (a[1]*b[2] - a[2]*b[1])
    end

    @inline function dot3(a, b)
        #=@inbounds=# return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
    end
end # module Rasterization
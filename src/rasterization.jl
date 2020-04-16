module Rasterization
    using Base.Threads
    using LinearAlgebra
    using Dates
    using Printf
    using NetCDF
    using Profile

    include("util.jl")

    struct TetrahedronFace
        n               :: NTuple{3, Float64}
        d               :: Float64
        dist_p_excl     :: Float64
    end

    mutable struct ThreadContext
        tet_points      :: Array{Float64, 2}
        tet_aabb        :: Array{Float64, 2}
        x_min           :: Int32
        y_min           :: Int32
        z_min           :: Int32
        x_max           :: Int32
        y_max           :: Int32
        z_max           :: Int32
        tet_dict        :: Array{Int32, 2}
        p               :: Array{Float64, 1}
        face_points     :: Array{SubArray, 1}
        n               :: Array{Float64, 1}
        tet_faces       :: Array{TetrahedronFace, 1}
        solved_problems :: Int
    end

    const EDGE_TOL = .0

    function rasterize(tetrahedra::AbstractArray{T, 2} where T <: Integer, 
                       points::AbstractArray{F, 2} where F <: AbstractFloat,
                       times::AbstractArray{Float64, 1}, 
                       sampling_rate::NTuple{3,Float64}, 
                       vars,
                       out_filename::AbstractString,
                       mem_limit::Int64=1024*1024*1024*8)

        domain_x = (minimum(points[1,:]), maximum(points[1,:]))
        domain_y = (minimum(points[2,:]), maximum(points[2,:]))
        domain_z = (minimum(points[3,:]), maximum(points[3,:]))

        println("Domain is $domain_x × $domain_y × $domain_z.")

        num_samples_x = ceil(Int32, (domain_x[end] - domain_x[1]) / sampling_rate[1])
        num_samples_y = ceil(Int32, (domain_y[end] - domain_y[1]) / sampling_rate[2])
        num_samples_z = ceil(Int32, (domain_z[end] - domain_z[1]) / sampling_rate[3])

        println("Sampling into $num_samples_x × $num_samples_y × $num_samples_z cells.")

        bytes_per_var_and_time = sizeof(Int64) * num_samples_x * num_samples_y
        bytes_per_grid = sizeof(UInt8) * num_samples_x * num_samples_y

        n_vars = 3
        n_times_per_iteration = (mem_limit - bytes_per_grid) ÷ (n_vars * bytes_per_var_and_time)
        n_iterations = ceil(Int32, length(times) / n_times_per_iteration)

        println("RAM limit: ", Util.human_readable_size(mem_limit), "; Can work on ", n_times_per_iteration, 
                " timesteps at once (therefore needing ", n_iterations , " iterations).")

        nccreate(out_filename, "u", 
                 "y", [domain_y[1] + i * sampling_rate[2] for i ∈ 1:num_samples_y], 
                 "x", [domain_x[1] + i * sampling_rate[1] for i ∈ 1:num_samples_x], 
                 "time", times, Dict("units"=>"seconds"))

        nccreate(out_filename, "v", 
                 "y", "x", "time")

        nccreate(out_filename, "w", 
                 "y", "x", "time")

        print("Allocating data structures and thread contexts... ")

        n_threads = nthreads()

        grid_locks = [ReentrantLock() for x ∈ 1:n_threads, y ∈ 1:n_threads]
        grid_running_avg_u = Array{Float64, 3}(undef, (num_samples_y, num_samples_x, n_times_per_iteration))
        grid_running_avg_v = Array{Float64, 3}(undef, (num_samples_y, num_samples_x, n_times_per_iteration))
        grid_running_avg_w = Array{Float64, 3}(undef, (num_samples_y, num_samples_x, n_times_per_iteration))
        grid_sample_cnts = Array{UInt16, 2}(undef, (num_samples_y, num_samples_x))

        print_lock = ReentrantLock()
        n_tetrahedra = size(tetrahedra, 2)

        total_problems = length(times) * n_tetrahedra

        println("Done.")

        start_time :: DateTime = now()
        last_printed_time = start_time
        print_interval = Second(2)
        n_tetrahedra_per_thread = n_tetrahedra / n_threads

        for iteration ∈ 1:n_iterations
            n_times = min(length(times) - (iteration - 1) * n_times_per_iteration, n_times_per_iteration)
            grid_running_avg_u[:,:,1:n_times] .= 0
            grid_running_avg_v[:,:,1:n_times] .= 0
            grid_running_avg_w[:,:,1:n_times] .= 0
            grid_sample_cnts .= 0

            Threads.@threads for thread_id ∈ 1:n_threads
                ctxt = ThreadContext(
                    Array{Float64, 2}(undef, (4, 3)),
                    Array{Float64, 2}(undef, (4, 2)),
                    0, 0, 0,
                    0, 0, 0,
                    Array{UInt8, 2}(undef, (64, 64)),
                    Array{Float64, 1}(undef, 3),
                    Array{SubArray, 1}(undef, 4),
                    Array{Float64, 1}(undef, 3),
                    Array{TetrahedronFace, 1}(undef, 4),
                    0)

                thread_range_min = floor(Int32, (thread_id - 1) * n_tetrahedra_per_thread) + 1
                thread_range_max = thread_id == n_threads ? n_tetrahedra : floor(Int32, thread_id * n_tetrahedra_per_thread)

                @profile begin
                for tet_id ∈ thread_range_min:thread_range_max
                    tetrahedron = (@view tetrahedra[:, tet_id])
                    for i ∈ 1:4, dim ∈ 1:3
                        ctxt.tet_points[i, dim] = points[dim, tetrahedron[i]]
                    end

                    for dim ∈ 1:3
                        ctxt.tet_aabb[dim, 1] = minimum(ctxt.tet_points[:, dim])
                        ctxt.tet_aabb[dim, 2] = maximum(ctxt.tet_points[:, dim])
                    end

                    ctxt.x_min = floor(Int32,(ctxt.tet_aabb[1,1] - domain_x[1]) / sampling_rate[1]) + 1
                    ctxt.y_min = floor(Int32,(ctxt.tet_aabb[2,1] - domain_y[1]) / sampling_rate[2]) + 1
                    ctxt.z_min = floor(Int32,(ctxt.tet_aabb[3,1] - domain_z[1]) / sampling_rate[3]) + 1

                    ctxt.x_max = ceil(Int32,(ctxt.tet_aabb[1,2] - domain_x[1]) / sampling_rate[1])
                    ctxt.y_max = ceil(Int32,(ctxt.tet_aabb[2,2] - domain_y[1]) / sampling_rate[2])
                    ctxt.z_max = ceil(Int32,(ctxt.tet_aabb[3,2] - domain_z[1]) / sampling_rate[3])

                    if thread_id == 1
                        t_now = now()
                        if (t_now - last_printed_time) >= print_interval && ctxt.solved_problems > 0
                            solved_problems = ctxt.solved_problems

                            last_printed_time = t_now
                            etr = (t_now-start_time) * (total_problems ÷ n_threads - solved_problems) ÷ (solved_problems)
                            hh = floor(etr, Dates.Hour)
                            etr -= floor(hh, Dates.Millisecond)
                            mm = floor(etr, Dates.Minute)
                            etr -= floor(mm, Dates.Millisecond)
                            ss = floor(etr, Dates.Second)

                            lock(print_lock)
                            @printf("Working on tetrahedron %d in timesteps %d-%d (%2.2f%% done). ETR: %02d:%02d:%02d    \r", 
                                          solved_problems, 
                                          (iteration - 1) * n_times_per_iteration + 1,
                                          (iteration - 1) * n_times_per_iteration + n_times,
                                          solved_problems/(total_problems/n_threads)*100,
                                          hh.value, mm.value, ss.value)
                            unlock(print_lock)
                        end
                    end

                    num_current_cells_x = ctxt.x_max - ctxt.x_min + 1
                    num_current_cells_y = ctxt.y_max - ctxt.y_min + 1
                    if size(ctxt.tet_dict, 1) < num_current_cells_y || size(ctxt.tet_dict, 2) < num_current_cells_x
                        ctxt.tet_dict = Array{UInt8, 2}(undef, (num_current_cells_y * 2, num_current_cells_x * 2))
                    end

                    for x ∈ 1:ctxt.x_max-ctxt.x_min+1, y ∈ 1:ctxt.y_max-ctxt.y_min+1
                        ctxt.tet_dict[y,x] = 0
                    end

                    for i ∈ 1:4
                        ctxt.face_points[i] = @view ctxt.tet_points[i,:]
                    end

                    function cross3!(a, b, ret)
                        ret[1] = (a[2]*b[3] - a[3]*b[2])
                        ret[2] = (a[3]*b[1] - a[1]*b[3])
                        ret[3] = (a[1]*b[2] - a[2]*b[1])
                    end

                    function dot3(a, b)
                        return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
                    end

                    # Array{(n, d, dist_p_excl)}
                    for i ∈ 1:4
                        # Exclude point i from the tetrahedron and consider the plane defined by the other 3 points
                        # If and only if (x,y,z) lies on the "positive" side (in normal direction) of all 4 planes, 
                        # it is inside the tetrahedron
                        p_excl =    ctxt.face_points[i]
                        p1 =        ctxt.face_points[(i%4)+1]
                        p2 =        ctxt.face_points[((i+1)%4)+1]
                        p3 =        ctxt.face_points[((i+2)%4)+1]
                    
                        # Calculate Hesse normal form of the plane defined by p1, p2, p3
                        # Normal vector
                        cross3!(p2-p1, p3-p1, ctxt.n)
                        normalize!(ctxt.n)
                    
                        # n1x1 + n2x2 + n3x3 + d = 0; solve for d
                        d = -dot3(p1, ctxt.n)
                        dist_p_excl = dot3(p_excl, ctxt.n) + d

                        # Ensure that p_excl is considered to be on the "positive" side of the plane
                        if dist_p_excl < 0
                            ctxt.n = -ctxt.n
                            d = -d
                            dist_p_excl = -dist_p_excl
                        end

                        ctxt.tet_faces[i] = TetrahedronFace((ctxt.n[1], ctxt.n[2], ctxt.n[3]), d, dist_p_excl)
                    end 

                    for x ∈ ctxt.x_min:ctxt.x_max
                        ctxt.p[1] = (x-1)*sampling_rate[1] + domain_x[1]
                        for y ∈ ctxt.y_min:ctxt.y_max
                            ctxt.p[2] = (y-1)*sampling_rate[2] + domain_y[1]
                            for z ∈ ctxt.z_min:ctxt.z_max
                                ctxt.p[3] = (z-1)*sampling_rate[3] + domain_z[1]
                            
                                accepted = true
                                # For each tetrahedron face, filter out cells that are not in the tetrahedron
                                for i ∈ 1:4
                                    face = ctxt.tet_faces[i]
                                    # p_excl is ALWAYS on the "positive" side of the plane. 
                                    # So, if p_excl and p lie on the same side, p is inside the tetrahedron
                                    dist_p = dot3(ctxt.p, face.n) + face.d
                                
                                    # The point lies approx. ON the face, reject in the tet on the lexicographically
                                    # more negative side of the face. This ensures that for two neighboring tetrahedra, the point is
                                    # only accepted for one of them.
                                    # Def.: lexicographical order: 
                                    #       a ≤ b ⟺ (a.x ≤ b.x) ∨ (a.x = b.x ∧ a.y ≤ b.y) ∨ (a.x = b.x ∧ a.y = b.y ∧ a.z ≤ b.z)
                                    # This is a total ordering on all vectors in R^3 and will therefore always be larger for the "positive"
                                    # normal vector in contrast to its negated counterpart.
                                    if (abs(dist_p) <= EDGE_TOL) 
                                        is_positive_side = 
                                            (ctxt.n[1] > 0.) || 
                                            (ctxt.n[1] == 0. && ctxt.n[2] > 0.) || 
                                            (ctxt.n[1] == 0. && ctxt.n[2] == 0. && ctxt.n[3] > 0.)

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
                                end

                                if accepted # accept cell and increment its 2D cell counter for the current tetrahedron by 1
                                    x_ = x - ctxt.x_min + 1; y_ = y - ctxt.y_min + 1
                                    ctxt.tet_dict[y_,x_] += 1
                                end
                            end
                        end
                    end


                    # Acquire grid region locks for writeback
                    lock_x_min = floor(Int32, ctxt.x_min / num_samples_x * n_threads) + 1
                    lock_x_max =  ceil(Int32, ctxt.x_max / num_samples_x * n_threads)
                    lock_y_min = floor(Int32, ctxt.y_min / num_samples_y * n_threads) + 1
                    lock_y_max =  ceil(Int32, ctxt.y_max / num_samples_y * n_threads)
                    for x ∈ lock_x_min:lock_x_max
                        for y ∈ lock_y_min:lock_y_max
                            lock(grid_locks[y, x])
                        end
                    end

                    for x ∈ 1:ctxt.x_max-ctxt.x_min+1, y ∈ 1:ctxt.y_max-ctxt.y_min+1
                        num_samples = ctxt.tet_dict[y, x]
                        if num_samples > 0
                            x_ = x + ctxt.x_min - 1; y_ = y + ctxt.y_min - 1
                            total_samples = (grid_sample_cnts[y_, x_] += num_samples)
                            for t ∈ 1:n_times
                                grid_running_avg_u[y_, x_, t] += (vars[t,1][tet_id]-grid_running_avg_u[y_, x_, t])*(num_samples/total_samples)
                                grid_running_avg_v[y_, x_, t] += (vars[t,2][tet_id]-grid_running_avg_v[y_, x_, t])*(num_samples/total_samples)
                                grid_running_avg_w[y_, x_, t] += (vars[t,3][tet_id]-grid_running_avg_w[y_, x_, t])*(num_samples/total_samples)
                            end
                        end
                    end

                    for x ∈ lock_x_max:-1:lock_x_min
                        for y ∈ lock_y_max:-1:lock_y_min
                            unlock(grid_locks[y, x])
                        end
                    end

                    ctxt.solved_problems += n_times
                end
                end#profile

                open("p-$thread_id.prof", "w") do p
                    Profile.print(p)
                end
            end

            start = [1, 1, (iteration - 1) * n_times_per_iteration + 1]
            count = [num_samples_y, num_samples_x, n_times]

            ncwrite(grid_running_avg_u, out_filename, "u", start=start, count=count)
            ncwrite(grid_running_avg_v, out_filename, "v", start=start, count=count)
            ncwrite(grid_running_avg_w, out_filename, "w", start=start, count=count)
        end

        println()
    end
end
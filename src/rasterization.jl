module Rasterization
    using Base.Threads
    using LinearAlgebra
    using Dates
    using Printf

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
    end

    function rasterize(tetrahedra, points, sampling_rate::NTuple{3,Float64})

        domain_x = (minimum(points[1,:]), maximum(points[1,:]))
        domain_y = (minimum(points[2,:]), maximum(points[2,:]))
        domain_z = (minimum(points[3,:]), maximum(points[3,:]))

        println("Domain is $domain_x × $domain_y × $domain_z.")

        num_samples_x = ceil(Int32, (domain_x[end] - domain_x[1]) / sampling_rate[1])
        num_samples_y = ceil(Int32, (domain_y[end] - domain_y[1]) / sampling_rate[2])
        num_samples_z = ceil(Int32, (domain_z[end] - domain_z[1]) / sampling_rate[3])

        println("Sampling into $num_samples_x × $num_samples_y × $num_samples_z cells.")

        print("Allocating data structures and thread contexts... ")

        grid_locks = [ReentrantLock() for x ∈ 1:8, y ∈ 1:8]
        grid = [Array{Tuple{Int32, UInt8}, 1}() for x ∈ 1:num_samples_x, y ∈ 1:num_samples_y]

        print_lock = ReentrantLock()
        incr_lock = ReentrantLock()
        avg_problem_size = [0., 0., 0.]
        solved_problems = 0
        n_threads = nthreads()
        n_tetrahedra = size(tetrahedra, 2)

        thread_contexts = [
            ThreadContext(
                Array{Float64, 2}(undef, (4, 3)),
                Array{Float64, 2}(undef, (4, 2)),
                0, 0, 0,
                0, 0, 0,
                Array{Int32, 2}(undef, (64, 64)),
                Array{Float64, 1}(undef, 3),
                Array{SubArray, 1}(undef, 4),
                Array{Float64, 1}(undef, 3),
                Array{TetrahedronFace, 1}(undef, 4)
            ) for i ∈ 1:n_threads
        ]

        println("Done.")

        start_time :: DateTime = now()
        print_counter = 0
        n_tetrahedra_per_thread = n_tetrahedra / n_threads

        Threads.@threads for thread_id ∈ 1:n_threads
            ctxt = thread_contexts[thread_id]

            lock(print_lock)
            thread_range_min = floor(Int32, (thread_id - 1) * n_tetrahedra_per_thread) + 1
            thread_range_max = thread_id == n_threads ? n_tetrahedra : floor(Int32, thread_id * n_tetrahedra_per_thread)
            @printf("Thread %2d working on tetrahedra {%d, ..., %d}\n", thread_id, thread_range_min, thread_range_max)
            unlock(print_lock)

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
                ctxt.z_min = floor(Int32,(ctxt.tet_aabb[3,1] - domain_z[1]) / sampling_rate[3])
                
                ctxt.x_max = ceil(Int32,(ctxt.tet_aabb[1,2] - domain_x[1]) / sampling_rate[1])
                ctxt.y_max = ceil(Int32,(ctxt.tet_aabb[2,2] - domain_y[1]) / sampling_rate[2])
                ctxt.z_max = ceil(Int32,(ctxt.tet_aabb[3,2] - domain_z[1]) / sampling_rate[3])
                
                if thread_id == 1
                    print_counter += 1
                    if print_counter % 2^11 == 0
                        lock(print_lock)
                        etr = (now()-start_time) * (n_tetrahedra - solved_problems) ÷ (solved_problems)
                        hh = floor(etr, Dates.Hour)
                        etr -= floor(hh, Dates.Millisecond)
                        mm = floor(etr, Dates.Minute)
                        etr -= floor(mm, Dates.Millisecond)
                        ss = floor(etr, Dates.Second)
                        @printf("Working on tetrahedron #%8d (%2.2f%% done). ETR: %02d:%02d:%02d\r", 
                                      solved_problems, 
                                      solved_problems/n_tetrahedra*100,
                                      hh.value, mm.value, ss.value)
                        unlock(print_lock)
                    end
                        solved_p1 = solved_problems+1
                        avg_problem_size[1] = (avg_problem_size[1] * (solved_problems/solved_p1) + (ctxt.x_max - ctxt.x_min + 1) / solved_p1)
                        avg_problem_size[2] = (avg_problem_size[2] * (solved_problems/solved_p1) + (ctxt.y_max - ctxt.y_min + 1) / solved_p1)
                        avg_problem_size[3] = (avg_problem_size[3] * (solved_problems/solved_p1) + (ctxt.z_max - ctxt.z_min + 1) / solved_p1)
                end
                
                num_current_cells_x = ctxt.x_max - ctxt.x_min + 1
                num_current_cells_y = ctxt.y_max - ctxt.y_min + 1
                if size(ctxt.tet_dict, 1) < num_current_cells_y || size(ctxt.tet_dict, 2) < num_current_cells_x
                    ctxt.tet_dict = Array{Int32, 2}(undef, (num_current_cells_y * 2, num_current_cells_x * 2))
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
                    ctxt.tet_faces[i] = TetrahedronFace((ctxt.n[1], ctxt.n[2], ctxt.n[3]), d, dist_p_excl)
                end 
                
                for x ∈ ctxt.x_min:ctxt.x_max
                    ctxt.p[1] = x*sampling_rate[1] + domain_x[1]
                    for y ∈ ctxt.y_min:ctxt.y_max
                        ctxt.p[2] = y*sampling_rate[1] + domain_y[1]
                        for z ∈ ctxt.z_min:ctxt.z_max
                            ctxt.p[3] = z*sampling_rate[1] + domain_z[1]
                        
                            # For each tetrahedron face, filter out cells that are not in the tetrahedron
                            for i ∈ 1:4
                                face = ctxt.tet_faces[i]
                                # p_excl is ALWAYS on the "positive" side of the plane. 
                                # So, if p_excl and p lie on the same side, p is inside the tetrahedron
                                dist_p = dot3(ctxt.p, face.n) + face.d
                            
                                # If signs of distance vectors do not match, p and p_excl are on different sides of the plane (reject)
                                if ((dist_p <= 0. && face.dist_p_excl > 0.) || (dist_p > 0. && face.dist_p_excl <= 0.))
                                    break
                                else # accept cell and increment its 2D cell counter for the current tetrahedron by 1
                                    x_ = x - ctxt.x_min + 1; y_ = y - ctxt.y_min + 1
                                    ctxt.tet_dict[y_,x_] += 1
                                end
                            end
                        end
                    end
                end
                

                # Acquire grid region locks for writeback
                lock_x_min = floor(Int32, ctxt.x_min / num_samples_x * 8) + 1
                lock_x_max =  ceil(Int32, ctxt.x_max / num_samples_x * 8)
                lock_y_min = floor(Int32, ctxt.y_min / num_samples_y * 8) + 1
                lock_y_max =  ceil(Int32, ctxt.y_max / num_samples_y * 8)
                for x ∈ lock_x_min:lock_x_max
                    for y ∈ lock_y_min:lock_y_max
                        lock(grid_locks[y, x])
                    end
                end

                i = 1
                for x ∈ 1:ctxt.x_max-ctxt.x_min+1, y ∈ 1:ctxt.y_max-ctxt.y_min+1
                    num_samples = ctxt.tet_dict[y, x]
                    if num_samples > 0
                        x_ = x + ctxt.x_min - 1; y_ = y + ctxt.y_min - 1
                        push!(grid[y_, x_], (tet_id, num_samples))
                        i += 1
                    end
                end

                for x ∈ lock_x_min:lock_x_max
                    for y ∈ lock_y_min:lock_y_max
                        unlock(grid_locks[y, x])
                    end
                end

                lock(incr_lock)
                solved_problems+=1
                unlock(incr_lock)
            end
        end

        println()
        println("Average problem size was: ", avg_problem_size, " cells per tetrahedron.")

        return grid, (num_samples_x, num_samples_y, num_samples_z)
    end
end
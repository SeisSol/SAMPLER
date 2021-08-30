module SAMPLER
    include("io/util.jl")
    include("io/args.jl")
    include("io/netcdf.jl")
    include("rasterizer/rasterization.jl")

    using NetCDF
    using Printf
    using Base.Filesystem
    using Base
    using SeisSolXDMF
    using Distributed, ClusterManagers
    using ..Args
    using ..Rasterization
    using ..Util
    using ..NC

    const WORKER_POOL = try
            n_workers = parse(Int, ENV["SLURM_NTASKS"])
            worker_procs = addprocs(SlurmManager(n_workers), exclusive=""; exeflags="--project")

            @everywhere include("io/util.jl")
            @everywhere include("io/args.jl")
            @everywhere include("rasterizer/rasterization.jl")
            @everywhere using ..Rasterization

            WorkerPool(worker_procs)
        catch _
            nothing
        end
    
    function main()
        println("Using $(Threads.nthreads()) threads per node.")

        ARGS = read_args()

        #============================================#
        # Read and preprocess args
        #============================================#

        mem_limit       = ARGS["memory-limit"]
        sampling_rate   = ARGS["sampling-rate"]
        domain          = ARGS["domain"]

        has_3d          = !isempty(ARGS["input-file-3d"])
        seafloor_only   = ARGS["seafloor-only"]
        
        seafloor_vars   = ARGS["seafloor-vars"]
        surface_vars    = ARGS["surface-vars"]
        volumetric_vars = has_3d ? ARGS["volumetric-vars"] : Dict{String,String}()
        
        water_height    = ARGS["water-height"]
        has_tanioka     = ARGS["tanioka"]

        is_parallel     = ARGS["parallel"]

        if is_parallel && isnothing(WORKER_POOL)
            error("Running in parallel but could not find SLURM job environment variables.")
        elseif is_parallel
            println("Running in distributed mode with $(length(WORKER_POOL)) workers.")
        else
            println("Running in serial mode.")
        end

        out_filename    = ARGS["output-file"]
        endswith(out_filename, ".nc") || (out_filename = out_filename * ".nc")
        println("Output will be written to \"$out_filename\".")

        surface_output  = !seafloor_only
        volume_output   = !seafloor_only && has_3d


        #============================================#
        # Compare timesteps of 2D and 3D files.
        # They must be equal.
        #============================================#

        xdmf2d = XDMFFile(ARGS["input-file-2d"])
        times  = timesteps_of(xdmf2d)

        if volume_output
            xdmf3d = XDMFFile(ARGS["input-file-3d"])
            times_3d = timesteps_of(xdmf3d)

            if (times_3d != times)
                throw(ArgumentError("Timesteps of 2D and 3D input files do not match!"))
            end

            xdmf3d = nothing
            times_3d = nothing
        end

        #============================================#
        # Compute timesteps to be processed from args
        #============================================#

        t_start = 1
        t_end   = length(times)

        if !isnothing(ARGS["output-time"])
            out_time = ARGS["output-time"]
            while t_start < length(times) && times[t_start + 1] ≤ out_time.t_start
                t_start += 1
            end

            while t_end > t_start && times[t_end - 1] ≥ out_time.t_end
                t_end -= 1
            end
        elseif !isnothing(ARGS["output-steps"])
            out_steps = ARGS["output-steps"]
            t_start   = min(1 + out_steps.t_start, length(times))
            t_end     = min(1 + out_steps.t_end, length(times))
        end

        #============================================#
        # Create NetCDF output file
        #============================================#

        empty_mapping    = Args.VarMapping()
        stat_var_mapping = Args.VarMapping()
        # "b" is the only static variable. If it is found in the mapping, it is separated into the stat_var_mapping.
        if haskey(seafloor_vars, "b"); 
            stat_var_mapping["b"] = seafloor_vars["b"]
            delete!(seafloor_vars, "b")
        end

        simplices_2d, points_2d = grid_of(xdmf2d)

        contexts = [(xdmf2d, simplices_2d, points_2d, seafloor_vars, stat_var_mapping, Rasterization.z_floor, has_tanioka)]

        if surface_output
            push!(contexts, (xdmf2d, simplices_2d, points_2d, surface_vars, empty_mapping, Rasterization.z_surface, false))
        end

        if volume_output
            xdmf3d = XDMFFile(ARGS["input-file-3d"])
            simplices_3d, points_3d = grid_of(xdmf3d)
            
            push!(contexts, (xdmf3d, simplices_3d, points_3d, volumetric_vars, empty_mapping, Rasterization.z_all, false))
        end

        contexts = [RasterizationContext(xdmf, simps, points, sampling_rate, mem_limit, 
                                         dyn_vars, stat_vars, z_range, tanioka, water_height, 
                                         t_start, t_end, times, out_filename, domain) 
                    for (xdmf, simps, points, dyn_vars, stat_vars, z_range, tanioka) ∈ contexts]

        ctx = contexts[1]
        create_netcdf(out_filename, 
                      [ctx.domain[Y][MIN] + i * ctx.sampling_rate[Y] for i ∈ 1:ctx.samples[Y]],
                      [ctx.domain[X][MIN] + i * ctx.sampling_rate[X] for i ∈ 1:ctx.samples[X]],
                      ctx.times[ctx.t_begin:ctx.t_end],
                      [stat_var_mapping], [seafloor_vars, surface_vars, volumetric_vars])

        rasterizer_params = [(ctx, bin(ctx)) for ctx ∈ contexts]

        #============================================#
        # Set up parallel workers (if enabled)
        #============================================#

        if is_parallel
            parallel_rasterizer_params = []

            get_t_start(ctx, iter) = ctx.t_begin + (iter - 1) * ctx.n_timesteps_per_iteration

            for (ctx, bins) ∈ rasterizer_params
                for iter ∈ 1:ctx.n_iterations
                    iter_t_start = get_t_start(ctx, iter)
                    iter_t_end = min(ctx.t_end, get_t_start(ctx, iter + 1) - 1)

                    iter_ctx = RasterizationContext(ctx.xdmf,
                        ctx.simplices, ctx.points, ctx.sampling_rate, mem_limit, 
                        ctx.dyn_var_mapping, ctx.stat_var_mapping, ctx.z_range, ctx.tanioka, ctx.water_height, 
                        iter_t_start, iter_t_end, times, ctx.out_filename, ctx.domain[X:Y]
                    )

                    push!(parallel_rasterizer_params, (iter_ctx, bins))
                end
            end

            rasterizer_params = nothing

            function rasterize_wrapper(ctx_bins_tuple)
                return rasterize(ctx_bins_tuple...)
            end

            pmap(rasterize_wrapper, WORKER_POOL, parallel_rasterizer_params)
        else
            for (ctx, bins) ∈ rasterizer_params
                rasterize(ctx, bins)
            end
        end
    end # main
end
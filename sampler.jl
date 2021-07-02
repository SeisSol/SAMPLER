#=
sampler.jl:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

include("io/util.jl")
include("io/args.jl")
include("io/xdmf.jl")
include("io/netcdf.jl")
include("rasterizer/rasterization.jl")

using Base.Threads
using Mmap
using NetCDF
using Profile
using Printf
using Base.Filesystem
using Pkg
using Base

println(Base.ARGS)
const ARGS = Args.read_args()

function main()

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

    out_filename    = ARGS["output-file"]
    endswith(out_filename, ".nc") || (out_filename = out_filename * ".nc")

    surface_output  = !seafloor_only
    volume_output   = !seafloor_only && has_3d

    #============================================#
    # Compare timesteps of 2D and 3D files.
    # They must be equal.
    #============================================#

    xdmf2d = XDMF.XDMFFile(ARGS["input-file-2d"])
    times  = XDMF.timesteps_of(xdmf2d)

    if has_3d
        xdmf3d = XDMF.XDMFFile(ARGS["input-file-3d"])
        times_3d = XDMF.timesteps_of(xdmf3d)

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
    # Process 2D seafloor
    #============================================#

    empty_mapping    = Args.VarMapping()
    stat_var_mapping = Args.VarMapping()
    # "b" is the only static variable. If it is found in the mapping, it is separated into the stat_var_mapping.
    if haskey(seafloor_vars, "b"); 
        stat_var_mapping["b"] = seafloor_vars["b"]
        delete!(seafloor_vars, "b")
    end

    Rasterization.rasterize(xdmf2d, seafloor_vars, stat_var_mapping, times, sampling_rate, out_filename, mem_limit;
                            create_nc_dyn_vars=[seafloor_vars, surface_vars, volumetric_vars],
                            create_nc_stat_vars=[stat_var_mapping],
                            z_range=Rasterization.z_floor, t_begin=t_start, t_end=t_end, water_height=water_height, 
                            tanioka=has_tanioka, domain=domain)

    GC.gc(true)

    #============================================#
    # Process 2D sea surface
    #============================================#

    if surface_output
        Rasterization.rasterize(xdmf2d, surface_vars, empty_mapping, times, sampling_rate, out_filename, mem_limit; 
                                z_range=Rasterization.z_surface, t_begin=t_start, t_end=t_end, water_height=water_height,
                                domain=domain)
    end


    xdmf2d = nothing
    GC.gc(true)

    #============================================#
    # Process 3D mesh
    #============================================#

    if volume_output
        xdmf3d = XDMF.XDMFFile(ARGS["input-file-3d"])

        Rasterization.rasterize(xdmf3d, volumetric_vars, empty_mapping, times, sampling_rate, out_filename, mem_limit; 
                                z_range=Rasterization.z_all, t_begin=t_start, t_end=t_end, water_height=water_height,
                                domain=domain)
    end
end

println("Using $(nthreads()) threads.")

Pkg.precompile()
println("Done.")

# Closest thing we have to __name__ == "__main__" in Julia :/
if !isinteractive()
    main()
end

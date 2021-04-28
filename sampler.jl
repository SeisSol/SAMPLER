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

const ARGS = Args.read_args()

function main()
    
    sampling_rate = ARGS["sampling-rate"]
    has_3d = !isempty(ARGS["input-file-3d"])
    has_tanioka = ARGS["tanioka"]
    seafloor_only = ARGS["seafloor-only"]

    seafloor_vars = ARGS["seafloor-vars"]
    surface_vars = ARGS["surface-vars"]
    volumetric_vars = ARGS["volumetric-vars"]

    #============================================#
    # Compare timesteps of 2D and 3D files.
    # They must be equal.
    #============================================#

    times    = XDMF.timesteps_of(ARGS["input-file-2d"])
    t_start  = ARGS["output-time"].t_start
    t_end    = ARGS["output-time"].t_end

    if has_3d
        times_3d = XDMF.timesteps_of(ARGS["input-file-3d"])

        if (times_3d != times)
            throw(ArgumentError("Timesteps of 2D and 3D input files do not match!"))
        end

        times_3d = nothing 
    end

    timestep_begin = 1
    while timestep_begin != length(times) && times[timestep_begin + 1] ≤ t_start
        timestep_begin += 1
    end

    timestep_end = length(times)
    while timestep_end != 1 && times[timestep_end - 1] ≥ t_end
        timestep_end -= 1
    end

    water_height = ARGS["water-height"]

    # Delete output file if it already exists
    out_filename = ARGS["output-file"]
    endswith(out_filename, ".nc") || (out_filename = out_filename * ".nc")

    surface_output = !seafloor_only
    volume_output = !seafloor_only && has_3d

    #============================================#
    # Process 2D seafloor
    #============================================#

    triangles,  points_2d = XDMF.grid_of(ARGS["input-file-2d"])

    # If tanioka, add horizontal seafloor displacements to seafloor_vars if not already in there
    if has_tanioka
        for var ∈ ["U", "V"]; if !haskey(seafloor_vars, var); seafloor_vars[var] = var; end; end
    end

    # Ensure bathymetry output
    if !haskey(seafloor_vars, "b"); seafloor_vars["b"] = "b"; end

    all_out_names = [collect(values(seafloor_vars)); 
                     surface_output ? collect(values(surface_vars)) : []; 
                     volume_output ? collect(values(volumetric_vars)) : []]
    
    # b is the only static variable, this list only covers dynamic ones
    all_out_names = filter(name -> name != seafloor_vars["b"], all_out_names)

    in_names = collect(keys(seafloor_vars))
    Rasterization.rasterize(triangles, points_2d, XDMF.data_of(ARGS["input-file-2d"], in_names...), in_names, seafloor_vars, 
                            times, sampling_rate, out_filename, ARGS["memory-limit"], 
                            z_range=Rasterization.z_floor, create_file_vars=all_out_names, 
                            t_begin=timestep_begin, t_end=timestep_end, 
                            water_height=water_height, tanioka=has_tanioka)

    GC.gc(true)

    #============================================#
    # Process 2D sea surface
    #============================================#

    if surface_output
        in_names = collect(keys(surface_vars))
        Rasterization.rasterize(triangles, points_2d, XDMF.data_of(ARGS["input-file-2d"], in_names...), in_names, surface_vars,
                                times, sampling_rate, out_filename, ARGS["memory-limit"], 
                                z_range=Rasterization.z_surface, 
                                t_begin=timestep_begin, t_end=timestep_end, water_height=water_height)
    end


    triangles = nothing
    points_2d = nothing
    GC.gc(true)

    #============================================#
    # Process 3D mesh
    #============================================#

    if volume_output
        tetrahedra, points_3d = XDMF.grid_of(ARGS["input-file-3d"])

        in_names = collect(keys(volumetric_vars))
        Rasterization.rasterize(tetrahedra, points_3d, XDMF.data_of(ARGS["input-file-3d"], in_names...), in_names, volumetric_vars,
                                times, sampling_rate, out_filename, ARGS["memory-limit"],
                                t_begin=timestep_begin, t_end=timestep_end, water_height=water_height)
    end
end

println("Using $(nthreads()) threads.")

Pkg.precompile()
println("Done.")

main()

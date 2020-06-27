#=
main:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

include("util.jl")
include("args.jl")
include("xdmf.jl")
include("kajiura.jl")
include("rasterization.jl")

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
    has_kajiura = ARGS["kajiura"]

    if has_3d && has_kajiura
        throw(ArgumentError("Cannot process 3D input file when Kajiura filter is active!"))
    end

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

    load_balancer = ARGS["load-balancer"]
    if load_balancer == "naive"
        load_balancer = Main.Rasterization.naive
    elseif load_balancer == "count"
        load_balancer = Main.Rasterization.count
    else
        load_balancer = Main.Rasterization.workload
    end

    # Delete output file if it already exists
    out_filename = ARGS["output-file"]
    endswith(out_filename, ".nc") || (out_filename = out_filename * ".nc")

    #============================================#
    # Process 2D seafloor
    #============================================#

    triangles,  points_2d = XDMF.grid_of(ARGS["input-file-2d"])

    Rasterization.rasterize(triangles, points_2d, XDMF.data_of(ARGS["input-file-2d"], "W"), ["W"], 
                            times, sampling_rate, out_filename, ARGS["memory-limit"], 
                            z_range=Rasterization.z_floor, create_file=true, kajiura=has_kajiura, 
                            t_begin=timestep_begin, t_end=timestep_end, load_balancer=load_balancer)

    GC.gc(true)

    #============================================#
    # Process 2D sea surface
    #============================================#

    if !has_kajiura
        Rasterization.rasterize(triangles, points_2d, XDMF.data_of(ARGS["input-file-2d"], "W"), ["W"], 
                                times, sampling_rate, out_filename, ARGS["memory-limit"], 
                                z_range=Rasterization.z_surface, 
                                t_begin=timestep_begin, t_end=timestep_end, load_balancer=load_balancer)
    end


    triangles = nothing
    points_2d = nothing
    GC.gc(true)

    #============================================#
    # Process 3D mesh
    #============================================#

    if has_3d
        tetrahedra, points_3d = XDMF.grid_of(ARGS["input-file-3d"])

        Rasterization.rasterize(tetrahedra, points_3d, XDMF.data_of(ARGS["input-file-3d"], "u", "v"), ["u", "v"], 
                                times, sampling_rate, out_filename, ARGS["memory-limit"], 
                                t_begin=timestep_begin, t_end=timestep_end, load_balancer=load_balancer)
    end
end

println("Using $(nthreads()) threads.")

Pkg.precompile()
println("Done.")

main()

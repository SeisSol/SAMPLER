#=
main:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

include("util.jl")
include("args.jl")
include("xdmf.jl")
include("rasterization.jl")

using Base.Threads
using Mmap
using NetCDF
using Profile
using Printf
using Base.Filesystem
using Pkg

Pkg.precompile()
println("Done.")


function main()
    println("Using $(nthreads()) threads.")
    args = Args.read_args()

    sampling_rate     = args["sampling-rate"]

    times = XDMF.timesteps_of(args["input-file-3d"])
    times_2d = XDMF.timesteps_of(args["input-file-2d"])
    t_start  = args["output-time"].t_start
    t_end    = args["output-time"].t_end

    if (times != times_2d)
        throw(ArgumentError("Timesteps of 2D and 3D input files do not match!"))
    end

    times_2d = nothing 

    out_filename = args["output-file"]
    endswith(out_filename, ".nc") || (out_filename = out_filename * ".nc")

    triangles,  points_2d = XDMF.grid_of(args["input-file-2d"])

    Rasterization.rasterize(triangles, points_2d, XDMF.data_of(args["input-file-2d"], "W"), ["W"], 
                            times, sampling_rate, out_filename, args["memory-limit"]; z_range=Rasterization.z_floor, create_file=true)

    Profile.clear()
    GC.gc(true)

    Rasterization.rasterize(triangles, points_2d, XDMF.data_of(args["input-file-2d"], "W"), ["W"], 
                            times, sampling_rate, out_filename, args["memory-limit"]; z_range=Rasterization.z_surface)

    Profile.clear()

    triangles = nothing
    points_2d = nothing
    GC.gc(true)
    tetrahedra, points_3d = XDMF.grid_of(args["input-file-3d"])

    Rasterization.rasterize(tetrahedra, points_3d, XDMF.data_of(args["input-file-3d"], "u", "v"), ["u", "v"], 
                            times, sampling_rate, out_filename, args["memory-limit"])

    Profile.clear()
end

main()

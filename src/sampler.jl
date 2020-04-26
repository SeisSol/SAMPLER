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
    var_names_3d      = args["vars-3d"]
    var_names_surface = args["vars-surface"]
    var_names_floor   = args["vars-floor"]

    tetrahedra, points_3d, times_3d = XDMF.grid_of(args["input-file-3d"])
    triangles,  points_2d, times_2d = XDMF.grid_of(args["input-file-2d"])

    if (times_3d != times_2d)
        throw(ArgumentError("Timesteps of 2D and 3D input files do not match!"))
    end

    filename_3d      = splitext(args["output-file"])[1] * ".nc"
    filename_surface = splitext(args["output-file"])[1] * "-surface.nc"
    filename_floor   = splitext(args["output-file"])[1] * "-floor.nc"

    for fn ∈ [filename_3d, filename_floor, filename_surface]
        isfile(fn) && rm(fn)
    end

    Rasterization.rasterize(tetrahedra, points_3d, XDMF.data_of(args["input-file-3d"], var_names_3d...), var_names_3d, 
                            times_3d, sampling_rate, filename_3d, args["memory-limit"])

    Rasterization.rasterize(triangles, points_2d, XDMF.data_of(args["input-file-2d"], var_names_surface...), var_names_surface, 
                            times_2d, sampling_rate, filename_surface, args["memory-limit"])

    Rasterization.rasterize(triangles, points_2d, XDMF.data_of(args["input-file-2d"], var_names_floor...), var_names_floor, 
                            times_2d, sampling_rate, filename_floor, args["memory-limit"])
end

main()

#=
main:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

include("args.jl")
include("xdmf.jl")
include("rasterization.jl")

using Base.Threads
using Mmap
using NetCDF
using Profile
using Printf
using Base.Filesystem


function main()
    args = Args.read_args()

    sampling_rate = (75., 150., 100.) # (dx, dy, dz) in meters

    println("Using $(nthreads()) threads.")
    tetrahedra, points, times = XDMF.grid_of(args["input-file"])
    var_names = ["u", "v", "w"]
    if Filesystem.isfile(args["output-file"]); Filesystem.rm(args["output-file"]) end
    Rasterization.rasterize(tetrahedra, points, times, sampling_rate, XDMF.data_of(args["input-file"], var_names...), args["output-file"])
end

main()

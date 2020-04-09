#=
main:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

include("args.jl")
include("xdmf.jl")

using Base.Threads


function main()
    args = Args.read_args()

    println("Using $(nthreads()) threads.")
    XDMF.XDMFFile(args["input-file"], (100., 100., 100.))
end

main()

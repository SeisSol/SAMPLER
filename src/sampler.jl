#=
main:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

include("args.jl")

function main()
    args = Args.read_args()
    println(args)
end

main()

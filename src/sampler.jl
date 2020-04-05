#=
main:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

include("args.jl")
include("quadtrees.jl")
include("xdmf.jl")

function main()
    #args = Args.read_args()

    region = Quadtrees.Quad(-3., -4., 6., 8.)
    points = [
        Quadtrees.Point(-2., 1., "A"),
        Quadtrees.Point(2., 1., "B"),
        Quadtrees.Point(-2., -1., "C"),
        Quadtrees.Point(2., -1., "D"),
        Quadtrees.Point(-3., 7., "E"),
        Quadtrees.Point(3., 6., "F")
    ]

    qt = Quadtrees.Quadtree(region, points)

    XDMF.XDMFFile("C:/Users/Max/Downloads/output2/trunc_out.xdmf", (7500, 7500))
end

main()

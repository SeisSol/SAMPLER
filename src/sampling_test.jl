include("util.jl")
include("kajiura.jl")
include("rasterization.jl")
using Main.Rasterization

simplices = [[1,4,6,7] [1,2,4,6] [1,3,4,7] [1,6,5,7] [8,7,6,4]]
points = [
    [0., 0., 0.] [1., 0., 0.] [0., 1., 0.] [1., 1., 0.] [0., 0., 1.] [1., 0., 1.] [0., 1., 1.] [1., 1., 1.]
    ]
times = [0.]
tet_vals = [1.,0.,0.,0.,0.]
in_vars = Array{AbstractArray, 2}(undef, (1, 2))
in_vars[1,1]=tet_vals
in_vars[1,2]=tet_vals

for δ ∈ 2:200
    Rasterization.rasterize(
        simplices,
        points,
        in_vars,
        ["u", "v"],
        [δ * 1. - 2.],
        (1. / δ,1. / δ,1. / δ),
        "sampling-test-$(δ).nc",
        1000000000,
        create_file=true
    )
end

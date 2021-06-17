include("../io/util.jl")
include("../io/args.jl")
include("../io/netcdf.jl")
include("../io/xdmf.jl")
include("../rasterizer/rasterization.jl")

using Test

ATOL=0.05

function norm(vec)
    mag = sqrt(sum(abs2, vec))
    return Tuple(el / mag for el ∈ vec)
end

function make_context(simplices, points; 
                      sampling_rate=(100., 100., 100.), 
                      mem_limit=4*1024^3, 
                      dyn_var_mapping=Dict("W"=>"d"),
                      stat_var_mapping=Dict("b"=>"b"),
                      z_range=Rasterization.z_all,
                      tanioka=false,
                      water_height=0.,
                      t_begin=0,
                      t_end=0,
                      times=[0.],
                      out_filename="temp.nc")
                      
    Rasterization.RasterizationContext(simplices, points, sampling_rate, mem_limit, 
                                       dyn_var_mapping, stat_var_mapping, z_range, tanioka, 
                                       water_height, t_begin, t_end, times, out_filename)
end

@testset "Rasterization" begin

    simplices2D = [[1, 2, 3   ] [2, 3, 4   ]]
    simplices3D = [[1, 2, 3, 4] [2, 3, 4, 5]]

    A = [1., 1., 1.]
    B = [2., 1., 2.]
    C = [1., 2., 1.]
    D = [2., 1., 1.]
    E = [2., 2., 2.]
    points2D = [A B C   E]
    points3D = [A B C D E]

    dimensions = Dict(2=>(simplices=simplices2D, points=points2D),
                     3=>(simplices=simplices3D, points=points3D))

    # Each case is [Dimension, Simplex Index, Simplex HNFs]
    hnf_test_cases = [
        (dim=2, 
         simp_idx=1, 
         hnfs=[Rasterization.HesseNormalForm(norm((-1., -1.,  0.)), -2.12),
               Rasterization.HesseNormalForm(     ( 1.,  0.,  0.) ,  1.00),
               Rasterization.HesseNormalForm(     ( 0.,  1.,  0.) ,  1.00)]),

        (dim=3, 
         simp_idx=1, 
         hnfs=[Rasterization.HesseNormalForm(norm((-1., -1.,  0.)), -2.12),
               Rasterization.HesseNormalForm(     ( 0.,  0.,  1.) ,  1.00),
               Rasterization.HesseNormalForm(     ( 0.,  1.,  0.) ,  1.00),
               Rasterization.HesseNormalForm(norm(( 1.,  0., -1.)),  0.00)])
    ]

    @testset "calculate_hnfs! ($(test_case.dim)D, simplex $(test_case.simp_idx))" for test_case ∈ hnf_test_cases
        dimension = dimensions[test_case.dim]
        ctx = make_context(dimension.simplices, dimension.points)
        l_simp_points = [@view ctx.points[:, p] for p ∈ ctx.simplices[:, test_case.simp_idx]]
        l_simp_hnfs = Array{Rasterization.HesseNormalForm, 1}(undef, ctx.n_simplex_points)

        Rasterization.calculate_hnfs!(ctx, l_simp_points, l_simp_hnfs)

        for p ∈ 1:ctx.n_simplex_points
            @test test_case.hnfs[p].d ≈ l_simp_hnfs[p].d atol=ATOL

            for dim ∈ 1:3
                @test test_case.hnfs[p].n0[dim] ≈ l_simp_hnfs[p].n0[dim] atol=ATOL
            end
        end
    end

    is_in_test_cases = [
        (dim=2, 
        simp_idx=1,
        points=[A, B, C, (A .+ B .+ C) / 3., A .- [1., .5, .2], B .+ [2., 0., 0.], C .+ [.1, .2, 0.]], 
        is_in=[true, false, false, true, false, false, false]) # B and C intersect an which has a negative normal. Thus, they are not in the triangle themselves

        (dim=3, 
        simp_idx=1,
        points=[A, B, C, D, (A .+ B .+ C .+ D) / 4., (A .+ C .+ D) / 3., A .- [1., .5, .2], B .- [.5, .5, 0.], D .+ [0., .2, 0.]], 
        is_in=[true, false, false, false, true, true, false, false, false])
    ]

    @testset "is_in_simplex ($(test_case.dim)D, simplex $(test_case.simp_idx))" for test_case ∈ is_in_test_cases
        dimension = dimensions[test_case.dim]
        ctx = make_context(dimension.simplices, dimension.points)
        l_simp_points = [@view ctx.points[:, p] for p ∈ ctx.simplices[:, test_case.simp_idx]]
        l_simp_hnfs = Array{Rasterization.HesseNormalForm, 1}(undef, ctx.n_simplex_points)

        Rasterization.calculate_hnfs!(ctx, l_simp_points, l_simp_hnfs)

        for (point, result) ∈ zip(test_case.points, test_case.is_in)
            @test result == Rasterization.is_in_simplex(ctx, l_simp_hnfs, point)
        end
    end
    
    @testset "is_bathy" begin end
    Rasterization.is_bathy
    @testset "bin" begin end
    Rasterization.bin
    @testset "iterate" begin end
    Rasterization.iterate
    @testset "iterate_simps!" begin end
    Rasterization.iterate_simps!
    @testset "cross3!" begin end
    Rasterization.cross3!
    @testset "dot3" begin end
    Rasterization.dot3
    @testset "bary_coords_2d" begin end
    Rasterization.bary_coords_2d
end


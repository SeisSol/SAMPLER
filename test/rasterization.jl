include("../src/io/util.jl")
include("../src/io/args.jl")
include("../src/io/netcdf.jl")
include("../src/rasterizer/rasterization.jl")

using Test
using Base.Threads

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
                      out_filename="temp.nc",
                      custom_domain=((-Inf, Inf), (-Inf, Inf)),
                      max_threads=typemax(Int))
                      
    Rasterization.RasterizationContext(nothing, simplices, points, sampling_rate, mem_limit, 
                                       dyn_var_mapping, stat_var_mapping, z_range, tanioka, 
                                       water_height, t_begin, t_end, times, out_filename, custom_domain, max_threads=max_threads)
end

@testset "Rasterization" begin
    if nthreads() == 1
        error("Binning tests only work in a multi-threaded environment. Please run with at least 3 threads.")
    end

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

    #============================================#
    # calculate_hnfs!
    #============================================#

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

    #============================================#
    # is_in_simplex
    #============================================#

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
    
    #============================================#
    # is_bathy
    #============================================#

    is_bathy_test_cases = [
        (bathy_heights=[1998., 2001., 2001.], water_height=2000., is_bathy=true),
        (bathy_heights=[-20., -19., -17.], water_height=0., is_bathy=true),
        (bathy_heights=[-20., -20., -20.], water_height=0., is_bathy=true),
        (bathy_heights=[0., 0., 0.], water_height=20., is_bathy=true),
        (bathy_heights=[0., 0., 0.], water_height=0., is_bathy=false),
        (bathy_heights=[2000., 2000., 2000.], water_height=2000., is_bathy=false),
        (bathy_heights=[2001., 2001., 2001.], water_height=2000., is_bathy=true)
    ]

    @testset "is_bathy (wh=$(test_case.water_height), b=$(test_case.bathy_heights))" for test_case ∈ is_bathy_test_cases
        bathys = test_case.bathy_heights
        points = [[1., 1., bathys[1]] [4., 4., bathys[2]] [2., 2., bathys[3]]]
        @test test_case.is_bathy == Rasterization.is_bathy(points, test_case.water_height)
    end

    #============================================#
    # bin
    #============================================#

    nBins = 3
    cX = nBins*3
    cY = 9
    cZ = 1
    points = zeros((3, (cX+1)*(cY+1)*(cZ+1)))
    for i ∈ 0:size(points,2)-1
        points[:, i+1] = [i % (cX+1), i ÷ (cX+1) % (cY+1), i ÷ ((cX+1) * (cY+1)) % (cZ+1)]
    end

    tup2idx = (x, y, z) -> (x + (cX+1) * (y + (cY+1) * z) + 1)

    for z ∈ 0:cZ, y ∈ 0:cY, x ∈ 0:cX
        p = points[:,tup2idx(x,y,z)]
        if tuple(p...) != (x,y,z)
            println("$((x,y,z)) != $p")
        end
    end

    simplices2D = []
    bin_ids = []
    z = 0
    water_height = 1.
    for z ∈ 0:cZ, y ∈ 0:0, x ∈ 0:cX-1
        # For each (x,y,z), which we define as the left upper point of the triangle, add two triangles:
        # (x,y  ,z)  o---o (x+1,y  ,z)
        #            |\ 1|
        #            | \ |
        #            |2 \|
        # (x,y+1,z)  o---o (x+1,y+1,z)
        push!(simplices2D, [tup2idx(x, y, z), tup2idx(x+1, y, z), tup2idx(x+1, y+1, z)])
        push!(simplices2D, [tup2idx(x, y, z), tup2idx(x+1, y+1, z), tup2idx(x, y+1, z)])
        bin_id = z < water_height ? 1 + x ÷ (cX ÷ nBins) : 2 * nBins + 1
        push!(bin_ids, bin_id)
        push!(bin_ids, bin_id)
    end

    for z ∈ 0:cZ, y ∈ 0:0, xBin ∈ 1:nBins-1
        # Triangles all span one bucket border
        x=(cX ÷ nBins)*xBin-1
        push!(simplices2D, [tup2idx(x, y, z), tup2idx(x+2, y, z), tup2idx(x+2, y+1, z)])
        push!(simplices2D, [tup2idx(x, y, z), tup2idx(x+2, y+1, z), tup2idx(x, y+1, z)])
        bin_id = z < water_height ? nBins + xBin : 2 * nBins + 1
        push!(bin_ids, bin_id)
        push!(bin_ids, bin_id)
    end

    for z ∈ 0:cZ, y ∈ 0:0
        # Triangles all span multiple borders
        push!(simplices2D, [tup2idx(1, y, z), tup2idx(cX, y, z), tup2idx(cX, y+1, z)])
        push!(simplices2D, [tup2idx(1, y, z), tup2idx(cX, y+1, z), tup2idx(1, y+1, z)])
        bin_id = z < water_height ? 2 * nBins : 2 * nBins + 1
        push!(bin_ids, bin_id)
        push!(bin_ids, bin_id)
    end

    simplices2D = [simplices2D[i][dim] for dim ∈ 1:3, i ∈ 1:length(simplices2D)]

    @testset "bin" begin
        ctx = make_context(simplices2D, points, sampling_rate=(1., 1., 1.), water_height=water_height, z_range=Rasterization.z_floor, max_threads=3)
        bins = Rasterization.bin(ctx)

        for i ∈ 1:size(simplices2D,2)
            @test bin_ids[i] == bins.bin_ids[i]
        end
    end

    @testset "iterate" begin end
    Rasterization.iterate
    @testset "iterate_simps!" begin end
    Rasterization.iterate_simps!
    @testset "cross3!" begin end
    Rasterization.cross3!
    @testset "dot3" begin end
    Rasterization.dot3

    #============================================#
    # bary_coords_2d
    #============================================#

    triangle = [[-2., 1., 14.] [1., 3., 13.] [-1., -1., 12.]]

    bary_test_cases = [
        (point=triangle[:, 1],                           bary=(  1.,   0.,   0.)),
        (point=triangle[:, 2],                           bary=(  0.,   1.,   0.)),
        (point=triangle[:, 3],                           bary=(  0.,   0.,   1.)),
        (point=sum(triangle,           dims=2)[:,1] / 3, bary=(1/3 , 1/3 , 1/3 )),
        (point=sum(triangle[:, 1:2],   dims=2)[:,1] / 2, bary=(1/2 , 1/2 ,   0.)),
        (point=sum(triangle[:, 2:3],   dims=2)[:,1] / 2, bary=(  0., 1/2 , 1/2 )),
        (point=sum(triangle[:, 1:2:3], dims=2)[:,1] / 2, bary=(1/2 ,   0., 1/2 )),

        (point=triangle[:, 1] + 2 * (triangle[:, 2] - triangle[:, 1]), bary=(-1.,  2., 0.)),
        (point=triangle[:, 1] + 2 * (triangle[:, 3] - triangle[:, 1]), bary=(-1.,  0., 2.)),
        (point=triangle[:, 2] + 2 * (triangle[:, 3] - triangle[:, 2]), bary=( 0., -1., 2.)),
    ]

    @testset "bary_coords_2d" for test_case ∈ bary_test_cases
        calculated_bary = Rasterization.bary_coords_2d(test_case.point, triangle)
        for dim ∈ 1:3
            @test test_case.bary[dim] ≈ calculated_bary[dim] atol=1e-6
        end
    end
end

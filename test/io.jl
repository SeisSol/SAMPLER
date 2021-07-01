using NetCDF: nc_close
include("../io/util.jl")
include("../io/args.jl")
include("../io/xdmf.jl")
include("../io/netcdf.jl")

using Test
using Main.Args
using Main.XDMF
using Main.NC
using NetCDF
using Base

@testset "IO" begin
    workdir = pwd()
    if endswith(workdir, "test")
        cd("..")
    end

    @testset "NetCDF" begin
        filename = "test_resources/test_output.nc"

        x_vals = -5:.5:5
        y_vals = -3:.25:3
        t_vals = 0:.1:2

        stat_maps = [Dict("b"=>"Bathy")]
        dyn_maps = [Dict("W"=>"BtmDisp", "w"=>"BtmVel"), Dict("W"=>"TopDisp", "u"=>"TopU", "v"=>"TopV", "w"=>"TopW")]

        nc = create_netcdf(filename, x_vals, y_vals, t_vals, stat_maps, dyn_maps)

        nc_vars = keys(nc.vars)
        expected_keys = []
        for mapping ∈ [stat_maps; dyn_maps]
            expected_keys = [expected_keys; collect(values(mapping))]
        end
        println(expected_keys)
        @test all(varname -> varname in nc_vars, expected_keys)
        @test nc["x"] == x_vals
        @test nc["y"] == y_vals
        @test nc["time"] == t_vals

        NetCDF.close(nc)
        nc = get_netcdf(filename)
        nc_vars = keys(nc.vars)
        @test all(varname -> varname in nc_vars, expected_keys)

        NetCDF.close(nc)
        rm(filename)
    end

    @testset "XDMF" begin
        filename = "test_resources/output-surface.xdmf"
        xdmf = XDMFFile(filename)
        timesteps = timesteps_of(xdmf)
        @test length(timesteps) == 201
        @test timesteps == [float(i) for i ∈ 0:200]

        simplices, points = grid_of(xdmf)
        @test size(simplices) == (3, 172616)
        @test size(points) == (3, 86468)

        simplices = nothing
        points = nothing

        vars = ["U", "V", "W", "u", "v", "w"]
        timeframe = (10, 15)
        var_dict = data_of(xdmf, timeframe[1], timeframe[2], vars)
        @test length(var_dict) == length(vars)
        @test length(var_dict[vars[1]]) == (timeframe[2] - timeframe[1] + 1) # both start and stop are inclusive
        @test maximum(var_dict[vars[1]][1]) > 0.
        @test minimum(var_dict[vars[1]][1]) < 0.
    end
end
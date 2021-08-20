using NetCDF: nc_close
include("../src/io/util.jl")
include("../src/io/args.jl")
include("../src/io/netcdf.jl")

using Test
using NetCDF
using Base
using Base.Filesystem

@testset "IO" begin
    workdir = pwd()
    if endswith(workdir, "test")
        cd("..")
    end

    @testset "NetCDF" begin
        filename = tempname(tempdir(); cleanup=true)

        x_vals = -5:.5:5
        y_vals = -3:.25:3
        t_vals = 0:.1:2

        stat_maps = [Dict("b"=>"Bathy")]
        dyn_maps = [Dict("W"=>"BtmDisp", "w"=>"BtmVel"), Dict("W"=>"TopDisp", "u"=>"TopU", "v"=>"TopV", "w"=>"TopW")]

        nc = NC.create_netcdf(filename, x_vals, y_vals, t_vals, stat_maps, dyn_maps)

        nc_vars = keys(nc.vars)
        expected_keys = []
        for mapping ∈ [stat_maps; dyn_maps]
            expected_keys = [expected_keys; collect(values(mapping))]
        end
        
        @test all(varname -> varname in nc_vars, expected_keys)
        @test nc["x"] == x_vals
        @test nc["y"] == y_vals
        @test nc["time"] == t_vals

        nc = nothing # Finalizer closes file
        nc = NC.get_netcdf(filename)
        nc_vars = keys(nc.vars)
        @test all(varname -> varname in nc_vars, expected_keys)
    end
end
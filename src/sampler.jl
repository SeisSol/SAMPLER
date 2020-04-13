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


function main()
    args = Args.read_args()

    sampling_rate = (100., 100., 100.) # (dx, dy, dz) in meters

    println("Using $(nthreads()) threads.")
    tetrahedra, points, times = XDMF.grid_of(args["input-file"])
    grid, (samples_x, samples_y, samples_z) = Rasterization.rasterize(tetrahedra, points, sampling_rate)

    tetrahedra = nothing
    points = nothing

    size_y, size_x = size(grid)
    cell_tets = Array{Tuple{Int16, Array{Tuple{Int32, Int16}}}}(undef, (size_y, size_x))

    nc_d_y    = NcDim("y", samples_y, values=[i * sampling_rate[2] for i ∈ 1:size_y]) # TODO: respect grid offset
    nc_d_x    = NcDim("x", samples_x, values=[i * sampling_rate[1] for i ∈ 1:size_x]) # TODO: respect grid offset
    nc_d_time = NcDim("time", length(times), values=times)

    var_names = ["u", "v", "w"]
    
    nc_dims = [nc_d_y, nc_d_x, nc_d_time]
    nc_vars = map(name -> NcVar(name, nc_dims, t=Float64), var_names)

    nc_file = NetCDF.create(args["output-file"], nc_vars)

    write_buf = Array{Float64, 2}(undef, (3, size_y))

    for timestep ∈ 1:1#length(times)
        vars = XDMF.data_of(args["input-file"], 1, var_names...)
        for x ∈ 1:size_x
            for y ∈ 1:size_y
                write_buf .= 0.
                grid_cell = grid[y, x]
                z_samples :: UInt8 = reduce((acc, tup) -> acc + tup[2], grid_cell, init=0)

                for (tet_id, tet_samples) ∈ grid_cell
                    for var_id ∈ 1:length(var_names)
                        write_buf[var_id, y] += (tet_samples / z_samples) * vars[var_id][tet_id]
                    end
                end

                empty!(grid_cell) # free the memory that is not needed anymore

                start = [1, x, timestep]
                count = [size_y, 1, 1]

                for var_id ∈ 1:length(var_names)
                    NetCDF.putvar(nc_vars[var_id], write_buf[var_id,:], start, count)
                end
            end

            @printf("Processed %4d of %4d grid columns (%2.2f%% done)\r", x, samples_x, x/samples_x*100)
        end

        NetCDF.sync(nc_file)
    end

    NetCDF.close(nc_file)
end

main()

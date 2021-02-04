#=
sampler.jl:
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

include("util.jl")
include("args.jl")
include("xdmf.jl")
include("netcdf.jl")
include("rasterization.jl")

using Base.Threads
using Mmap
using NetCDF
using Profile
using Printf
using Base.Filesystem
using Pkg

const ARGS = Args.read_args()

function main()
    
    sampling_rate = ARGS["sampling-rate"]
    has_3d = !isempty(ARGS["input-file-3d"])
    has_kajiura = ARGS["kajiura"]
    has_tanioka = ARGS["tanioka"]

    #============================================#
    # Compare timesteps of 2D and 3D files.
    # They must be equal.
    #============================================#

    times    = XDMF.timesteps_of(ARGS["input-file-2d"])
    t_start  = ARGS["output-time"].t_start
    t_end    = ARGS["output-time"].t_end

    if has_3d
        times_3d = XDMF.timesteps_of(ARGS["input-file-3d"])

        if (times_3d != times)
            throw(ArgumentError("Timesteps of 2D and 3D input files do not match!"))
        end

        times_3d = nothing 
    end

    timestep_begin = 1
    while timestep_begin != length(times) && times[timestep_begin + 1] ≤ t_start
        timestep_begin += 1
    end

    timestep_end = length(times)
    while timestep_end != 1 && times[timestep_end - 1] ≥ t_end
        timestep_end -= 1
    end

    #============================================#
    # Load balancing parameters.
    # The naive balancer yields the best runtimes.
    # Do not use the other ones.
    #============================================#

    lb_autotune = ARGS["lb-autotune"]

    load_balancer = ARGS["load-balancer"]
    lb_params = (0., 0.)

    if lb_autotune
        if !has_3d
            println("The workload-LB currently can only be autotuned on tetrahedra. Please supply a 3D mesh also.")
            exit(1)
        end

        println("Performing LB-autotune.")

        load_balancer = Main.Rasterization.naive
    elseif load_balancer == "naive"
        load_balancer = Main.Rasterization.naive
    elseif load_balancer == "count"
        load_balancer = Main.Rasterization.count
    else
        load_balancer = Main.Rasterization.workload

        if isfile("sampler_lb_params.txt")
            str_params = readline("sampler_lb_params.txt")
            ls_params = split(str_params, ';')
            lb_params = (parse(Float64, ls_params[1]), parse(Float64, ls_params[2]))
        else
            println("To use the workload-LB, please run first with your scenario settings and --lb-autotune")
            println("Using fallback (naive) LB instead!")
            load_balancer = Main.rasterization.naive
        end
    end

    water_height = ARGS["water-height"]

    # Delete output file if it already exists
    out_filename = ARGS["output-file"]
    endswith(out_filename, ".nc") || (out_filename = out_filename * ".nc")

    if !lb_autotune
        #============================================#
        # Process 2D seafloor
        #============================================#

        triangles,  points_2d = XDMF.grid_of(ARGS["input-file-2d"])

        Rasterization.rasterize(triangles, points_2d, XDMF.data_of(ARGS["input-file-2d"], "W"), ["W"], 
                                times, sampling_rate, out_filename, ARGS["memory-limit"], 
                                z_range=Rasterization.z_floor, create_file=true, 
                                t_begin=timestep_begin, t_end=timestep_end, load_balancer=load_balancer, 
                                lb_params=lb_params, water_height=water_height)

        GC.gc(true)

        #============================================#
        # Process 2D sea surface
        #============================================#

        if !has_kajiura
            var_names = has_tanioka ? ["W", "U", "V"] : ["W"]
            Rasterization.rasterize(triangles, points_2d, XDMF.data_of(ARGS["input-file-2d"], var_names...), var_names, 
                                    times, sampling_rate, out_filename, ARGS["memory-limit"], 
                                    z_range=Rasterization.z_surface, 
                                    t_begin=timestep_begin, t_end=timestep_end, load_balancer=load_balancer, 
                                    lb_params=lb_params, water_height=water_height, tanioka=has_tanioka)
        end


        triangles = nothing
        points_2d = nothing
        GC.gc(true)
    end

    #============================================#
    # Process 3D mesh
    #============================================#

    if has_3d && !has_kajiura
        tetrahedra, points_3d = XDMF.grid_of(ARGS["input-file-3d"])

        Rasterization.rasterize(tetrahedra, points_3d, XDMF.data_of(ARGS["input-file-3d"], "u", "v"), ["u", "v"], 
                                times, sampling_rate, out_filename, ARGS["memory-limit"], create_file=lb_autotune,
                                t_begin=timestep_begin, t_end=timestep_end, load_balancer=load_balancer, 
                                lb_params=lb_params, lb_autotune=lb_autotune, water_height=water_height)
    end
end

println("Using $(nthreads()) threads.")

Pkg.precompile()
println("Done.")

main()

#=
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

module Args
    using ArgParse
    using Main.Util

    export read_args
    export Timespan

    struct Timespan
        t_start :: Float64
        t_end   :: Float64
    end

    function ArgParse.parse_item(::Type{Timespan}, x::AbstractString)
        x = strip(x)

        if x == "all"
            return Timespan(0., Inf64)
        end

        if occursin(',', x)
            timestamps = split(x, ',')

            if length(timestamps) != 2
                throw(ArgumentError(""""$x" does not match the format "start,end"!"""))
            end

            # If empty, assume t_start to be 0 and t_end to be infinity
            t_start :: Float64 = 0
            t_end   :: Float64 = Inf64

            try
                timestamps[1] = strip(timestamps[1])
                if length(timestamps[1]) != 0
                    t_start = parse(Float64, strip(timestamps[1]))
                end
            catch
                throw(ArgumentError(""""$(timestamps[1])" is not a floating point value!"""))
            end

            try
                timestamps[2] = strip(timestamps[2])
                if length(timestamps[2]) != 0
                    t_end = parse(Float64, timestamps[2])
                end
            catch
                throw(ArgumentError(""""$(timestamps[2])" is not a floating point value!"""))
            end

            if (t_start < 0)
                    throw(ArgumentError("Timestamps have to be non-negative!"))
            end

            if t_start > t_end
                throw(ArgumentError("The start timestamp must be less than the end timestamp!"))
            end

            return Timespan(t_start, t_end)
        else
            try
                t = parse(Float64, x)

                if t < 0
                    throw(ArgumentError(""""$t" must be non-negative!"""))
                end

                return Timespan(t, t)
            catch
                throw(ArgumentError(""""$t" is not a floating point value!"""))
            end
        end
    end

    function ArgParse.parse_item(::Type{NTuple{3, Float64}}, x::AbstractString)
        components = split(x, ',')

        if (length(components) ∉ [1, 3])
            throw(ArgumentError(""""$x" has $(length(components)) components but can only have 1 or 3!"""))
        end

        try
            return tuple(map(x -> parse(Float64, x), components)...)
        catch
            throw(ArgumentError("""The components of "$x" cannot be parsed as numbers!"""))
        end
    end

    function read_args()
        parser_settings = ArgParseSettings()
        @add_arg_table! parser_settings begin
            "--output-file", "-o"
                help = "The filename (without extension) to use for the output."
                default = ""
            "--output-time", "-t"
                help = """The time range of the output. Format: "[start,]end" or "all".\n
                                \t start and end both have to follow the "[[hh:]mm:]ss" format.
                                The end timestamp is bounded by the last input timestamp.\n
                                \t By default, every timestep will be output."""
                arg_type = Timespan
                default = Timespan(-Inf64, Inf64)
            "--sampling-rate", "-r"
                help = """The size in meters of one cell edge in the resampled output.\n
                                \t Format: dx[,dy,dz]. If only dx is given, dy and dz will be set equal to dx.\n
                                \t Examples: 100; 50,70,80"""
                arg_type = NTuple{3, Float64}
                default = (100., 100., 100.)
            "--water-height"
                help = """Some bathymetry grids place the seafloor at z = 0. Set the water height such that it matches the sea surface's z-coordinate."""
                arg_type = Float64
                default = 0.
            "--memory-limit", "-m"
                help = """The maximum amount of RAM the script should use. This is only a SOFT limit!\n
                                \t Examples: 8G; 2T; 512M\n
                                \t IEC-Prefixes are used: 1K = 1KiB = 1024B, ..."""
                default = "8G"
            "--load-balancer"
                help = """The static load balancer to choose when assigning simplices to threads. Options are:\n
                                \t naive: Every thread works on a bin that is |D| / n_threads in size\n
                                \t [DO NOT USE] count: Every thread works on a bin that approximately contains n_simplices / n_threads simplices\n
                                \t [DO NOT USE] workload: Every thread has roughly the same ESTIMATED (!) workload (read: processing time) in its bin. Run with --lb-autotune first.
                                """
                default="naive"
            "--lb-autotune"
                help = """[DO NOT USE] Perform a test run with the naive load balancer and collect data to autotune the workload-LB. 
                    The found parameters will be saved automatically and the workload-LB can then be used."""
                action = :store_true
            "--kajiura"
                help = """Only rasterize the seafloor displacements. The outputs can be processed by kajiura.jl."""
                action = :store_true
            "input-file-2d"
                help = "The SeisSol 2D output in XDMF format. Use the output file with the _surface suffix."
                required = true
            "input-file-3d"
                help = "The SeisSol 3D output in XDMF format. Use the output file without the _surface suffix."
                default = ""
        end

        args = parse_args(parser_settings)
        validate_args!(args)

        return args
    end

    function validate_args!(args::Dict)
        args["memory-limit"] = Main.Util.parse_size(args["memory-limit"])

        if !(args["load-balancer"] in ["naive", "count", "workload"])
            throw(ArgumentError("Load Balancer '$(args["load-balancer"])' does not exist!"))
        end
        
        return args
    end
end

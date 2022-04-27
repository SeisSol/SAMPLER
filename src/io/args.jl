module Args
    using ArgParse
    using ..Util

    export read_args, Timespan, VarMapping, SamplingTuple, DomainSize


    VarMapping      = Dict{String, String}
    SamplingTuple   = NTuple{3, Float64}
    DomainSize      = NTuple{2, NTuple{2, Float64}} # ((xmin,xmax), (ymin,ymax))

    struct Timespan{T <: Number}
        t_start :: T
        t_end   :: T
    end

    function ArgParse.parse_item(::Type{Timespan{T}}, x::AbstractString) where T <: Number
        x = strip(x)

        t_start :: T = zero(T)
        t_end   :: T = typemax(T)

        if x == "all"
            return Timespan{T}(t_start, t_end)
        end

        if occursin(',', x)
            timestamps = split(x, ',')

            if length(timestamps) != 2
                throw(ArgumentError(""""$x" does not match the format "start,end"!"""))
            end

            try
                timestamps[1] = strip(timestamps[1])
                if length(timestamps[1]) != 0
                    t_start = parse(T, strip(timestamps[1]))
                end
            catch
                throw(ArgumentError(""""$(timestamps[1])" is not a floating point value!"""))
            end

            try
                timestamps[2] = strip(timestamps[2])
                if length(timestamps[2]) != 0
                    t_end = parse(T, timestamps[2])
                end
            catch
                throw(ArgumentError(""""$(timestamps[2])" is not a floating point value!"""))
            end

            if (t_start < zero(T))
                throw(ArgumentError("Timestamps have to be non-negative!"))
            end

            if t_start > t_end
                throw(ArgumentError("The start timestamp must be less than the end timestamp!"))
            end

            return Timespan{T}(t_start, t_end)
        else
            try
                t = parse(T, x)

                if t < zero(T)
                    throw(ArgumentError(""""$t" must be non-negative!"""))
                end

                return Timespan{T}(t, t)
            catch
                throw(ArgumentError(""""$t" is not a floating point value!"""))
            end
        end
    end

    function ArgParse.parse_item(::Type{VarMapping}, x::AbstractString)
        entries = Set(split(x, ','))
        non_mapping_entries = filter(entry -> !occursin("=>", entry), entries)
        mapping_entries = setdiff(entries, non_mapping_entries)
        mapping_entries = map(str -> split(str, "=>"), collect(mapping_entries))

        entries = merge(Dict(map(str_list   -> String(str_list[1])  => String(str_list[2]), mapping_entries)),
                        Dict(map(entry      -> String(entry)        => String(entry),       collect(non_mapping_entries))))
        return entries
    end

    function ArgParse.parse_item(::Type{SamplingTuple}, x::AbstractString)
        components = split(x, ',')

        if (length(components) ∉ [1, 3])
            throw(ArgumentError(""""$x" has $(length(components)) components but can only have 1 or 3!"""))
        end

        if length(components) == 1
            components = [components[1], components[1], components[1]]
        end

        try
            return tuple(map(x -> parse(Float64, x), components)...)
        catch
            throw(ArgumentError("""The components of "$x" cannot be parsed as numbers!"""))
        end
    end

    function ArgParse.parse_item(::Type{DomainSize}, x::AbstractString)
        floats = split(x, ',')

        if length(floats) != 4
            throw(ArgumentError(""""$x" is not in the format "xmin,xmax,ymin,ymax"!"""))
        end

        try
            floats = tuple(map(x -> parse(Float64, x), floats)...)
            return ((floats[1], floats[2]), (floats[3], floats[4]))
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
                                \t start and end are interpreted as timeSTAMPs, not timeSTEP indices!
                                   The end timestep is bounded by the last input timestep.\n
                                \t If only end is given, only that timestep will be rasterized.\n
                                \t If start or end is an empty string, it will be assumed to be 0 or infinity respectively.\n
                                \t By default, every timestep will be output."""
                arg_type = Timespan{Float64}
            "--output-steps", "-s"
                help = """The time range of the output. Format: "[start,]end" or "all".\n
                                \t start and end are interpreted as timeSTEP indices, not timeSTAMPs!
                                   The end timestep is bounded by the last input timestep.\n
                                \t If only end is given, only that timestep will be rasterized.\n
                                \t If start or end is an empty string, it will be assumed to be 0 or infinity respectively.\n
                                \t By default, every timestep will be output."""
                arg_type = Timespan{UInt64}
            "--domain", "-d"
                help = """Specify a custom domain which should be rasterized. Format: "xmin,xmax,ymin,ymax".\n
                                \t By default, the whole domain will be rasterized. The specified domain will be clipped to the actual domain."""
                arg_type = DomainSize
                default = ((-Inf, Inf), (-Inf, Inf))
            "--sampling-rate", "-r"
                help = """The size in meters of one cell edge in the resampled output.\n
                                \t Format: dx[,dy,dz]. If only dx is given, dy and dz will be set equal to dx.\n
                                \t Examples: 100; 50,70,80"""
                arg_type = SamplingTuple
                default = (100., 100., 100.)
            "--water-height"
                help = """Some bathymetry grids place the seafloor at z = 0. Set the water height such that it matches the sea surface's z-coordinate."""
                arg_type = Float64
                default = 0.
            "--memory-limit", "-m"
                help = """The maximum amount of RAM the script should use. This is only a SOFT limit!\n
                                \t [CAUTION] Currently, this limit is quite imprecise, do not specify more than ~half of the available RAM here!\n
                                \t Examples: 8G; 2T; 512M\n
                                \t IEC-Prefixes are used: 1K = 1KiB = 1024B, ..."""
                default = "8G"
            "--tanioka"
                help = """Apply Tanioka's method for converting horizontal displacements to vertical ones during rasterization.\n
                                \t If the variables u1 or u2 are not found in --seafloor-vars, they will be added automatically."""
                action = :store_true
            "--seafloor-vars"
                help = """The names of all variables that should be extracted from the 2D grid at the seafloor, separated by commas.\n
                                \t Note that the bathymetry "b" will always be extracted from geometry. Use "b=>NewName" to rename
                                   the output variable.\n
                                \t [CAUTION] You cannot have the same output name for multiple variables, 
                                   even across --seafloor-vars, --surface-vars and --volumetric-vars!\n
                                \t Format: [M1[,M2[,...]] where Mi is either a variable name or a mapping like "u3=>d".\n
                                \t Examples: u1,u2,u3; u1,u2,u3=>eta; u1=>X,u2=>Y"""
                arg_type = VarMapping
                default=Dict("u3"=>"d", "b"=>"b")
            "--surface-vars"
                help = """The names of all variables that should be extracted from the 2D grid at the sea surface, separated by commas.\n
                                \t The same formatting rules apply as for --seafloor-vars."""
                arg_type = VarMapping
                default=Dict("u3"=>"eta")
            "--volumetric-vars"
                help = """The names of all variables that should be extracted from the 3D grid, separated by commas.\n
                                \t The same formatting rules apply as for --seafloor-vars."""
                arg_type = VarMapping
                default=Dict("v1"=>"v1", "v2"=>"v2")
            "--seafloor-only"
                help = """Only rasterize the seafloor displacements. The outputs can be processed by kajiura.jl."""
                action = :store_true
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
        args["memory-limit"] = parse_size(args["memory-limit"])

        if !isnothing(args["output-time"]) && !isnothing(args["output-steps"])
           throw(ArgumentError("Setting both '--output-time' and '--output-steps' simultaneously is not allowed!"))
        end

        if !(args["load-balancer"] in ["naive", "count", "workload"])
            throw(ArgumentError("Load Balancer '$(args["load-balancer"])' does not exist!"))
        end
        
        return args
    end
end

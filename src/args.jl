#=
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-03-30
=#

module Args
    using ArgParse

    export read_args
    export Timespan

    struct Timespan
        t_start :: Float32
        t_end   :: Float32
    end

    function ArgParse.parse_item(::Type{Timespan}, x::AbstractString)
        x = strip(x)

        if x == "all"
            return Timespan(0., Inf32)
        end

        if occursin(',', x)
            timestamps = split(x, ',')

            if length(timestamps) != 2
                throw(ArgumentError(""""$x" does not match the format "start,end"!"""))
            end

            # If empty, assume t_start to be 0 and t_end to be infinity
            t_start :: Float32 = 0
            t_end   :: Float32 = Inf32

            try
                timestamps[1] = strip(timestamps[1])
                if length(timestamps[1]) != 0
                    t_start = parse(Float32, strip(timestamps[1]))
                end
            catch
                throw(ArgumentError(""""$(timestamps[1])" is not a floating point value!"""))
            end

            try
                timestamps[2] = strip(timestamps[2])
                if length(timestamps[2]) != 0
                    t_end = parse(Float32, timestamps[2])
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
                t = parse(Float32, x)

                if t < 0
                    throw(ArgumentError(""""$t" must be non-negative!"""))
                end

                return Timespan(t, t)
            catch
                throw(ArgumentError(""""$t" is not a floating point value!"""))
            end
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
                                \tstart and end both have to follow the "[[hh:]mm:]ss" format.
                                The end timestamp is bounded by the last input timestamp.\n
                                \tBy default, every timestep will be output."""
                arg_type = Timespan
                default = Timespan(-Inf64, Inf64)
            "input-file"
                help = "The SeisSol output in XDMF format. Use either the 3D or 2D (surface) output."
                required = true
        end

        args = parse_args(parser_settings)
        validate_args!(args)

        return args
    end

    function validate_args!(args::Dict)
        # TODO
        return args
    end
end

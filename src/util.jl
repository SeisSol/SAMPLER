module Util
    using NetCDF

    export human_readable_size, parse_size

    const ORDERS_OF_MAGNITUDE = "KMGTPEZY"

    function human_readable_size(bytes:: Integer) :: String
        order_of_magnitude = floor(Int, log(1024, bytes)) # Order of magnitude of bytes in base 1024 (0 for B, 1 for KiB, ...)
        order_of_magnitude = min(order_of_magnitude, length(ORDERS_OF_MAGNITUDE))
        bytes ÷= 1024^order_of_magnitude

        prefix = if order_of_magnitude == 0; ""
                 else ORDERS_OF_MAGNITUDE[order_of_magnitude] * "i" 
                 end

        return "$bytes $(prefix)B"
    end

    function parse_size(size:: String) :: Integer
        if isnothing(match(r"^[0-9]+[KMGTPEZY]$", size))
            throw(ArgumentError("The given size string '$size' is invalid. Valid examples: 8G, 512M"))
        end

        magnitude_letter = size[end]
        magnitude = parse(Int, size[1:end-1])
        order_of_magnitude = findfirst(magnitude_letter, ORDERS_OF_MAGNITUDE)

        return magnitude * 1024^order_of_magnitude
    end

    function get_or_create_netcdf(filename :: AbstractString; 
                                  create = false, x_vals = [], y_vals = [], t_vals = []) :: NcFile
        
        nc :: NcFile =
            if create
                isfile(filename) && rm(filename)

                if isempty(x_vals) || isempty(y_vals) || isempty(t_vals) 
                    throw(ArgumentError("Dimension values have to be passed when creating a NetCDF file!"))
                end

                x_atts = Dict("units" => "m")
                y_atts = Dict("units" => "m")
                t_atts = Dict("units" => "seconds")
                uv_atts = Dict("units" => "m/s")
                bdh_atts = Dict("units" => "m")

                x_dim = NcDim("x", x_vals, x_atts)
                y_dim = NcDim("y", y_vals, y_atts)
                t_dim = NcDim("time", t_vals, t_atts)

                b_var = NcVar("b", [y_dim, x_dim], atts=bdh_atts)

                u_var = NcVar("u", [y_dim, x_dim, t_dim], atts=uv_atts)
                v_var = NcVar("v", [y_dim, x_dim, t_dim], atts=uv_atts)
                d_var = NcVar("d", [y_dim, x_dim, t_dim], atts=bdh_atts)
                h_var = NcVar("h", [y_dim, x_dim, t_dim], atts=bdh_atts)

                NetCDF.create(filename, b_var, u_var, v_var, d_var, h_var)
            else
                NetCDF.open(filename, mode=NC_WRITE)
            end

        NetCDF.nc_set_fill(nc.ncid, NetCDF.NC_NOFILL, Ptr{Int32}(0))
        return nc
    end
end
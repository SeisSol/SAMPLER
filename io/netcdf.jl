#=
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-07-13
=#

module NC
    using NetCDF
    using Main.Args

    export get_or_create_netcdf

    """
    Create the NetCDF file `filename` with coordinate values `x_vals`, `y_vals` and `t_vals` or open an existing one.

    The format of the NetCDF file is chosen to exactly match the requirements of this program.
    """
    function get_or_create_netcdf(filename :: AbstractString; 
                                  create = false, x_vals = [], y_vals = [], t_vals = [],
                                  create_static_vars = ["b"],
                                  create_dynamic_vars = ["u", "v", "d", "eta"]) :: NcFile
        
        nc :: NcFile =
            if create
                isfile(filename) && rm(filename)

                if isempty(x_vals) || isempty(y_vals) || isempty(t_vals) 
                    throw(ArgumentError("Dimension values have to be passed when creating a NetCDF file!"))
                end

                t_atts = Dict("units" => "seconds")
                vel_atts = Dict("units" => "m/s")
                dist_atts = Dict("units" => "m")
                no_attrs = Dict()

                notfound_units = Dict(map(var_name -> no_attrs, [create_static_vars; create_dynamic_vars]))

                units = Dict(
                    "b"=>dist_atts,
                    "d"=>dist_atts,
                    "eta"=>dist_atts,
                    "u"=>vel_atts,
                    "v"=>vel_atts
                )

                units = merge(notfound_units, units)

                x_dim = NcDim("x", x_vals, dist_atts)
                y_dim = NcDim("y", y_vals, dist_atts)
                t_dim = NcDim("time", t_vals, t_atts)

                static_vars = map(name -> NcVar(name, [x_dim, y_dim], atts=units[name]), create_static_vars)
                dynamic_vars = map(name -> NcVar(name, [x_dim, y_dim, t_dim], atts=units[name]), create_static_vars)

                NetCDF.create(filename, static_vars..., dynamic_vars...)
            else
                NetCDF.open(filename, mode=NC_WRITE)
            end

        NetCDF.nc_set_fill(nc.ncid, NetCDF.NC_NOFILL, Ptr{Int32}(0))
        return nc
    end
end
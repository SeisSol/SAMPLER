#=
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-07-13
=#

module NC
    using NetCDF

    export get_or_create_netcdf

    """
    Create the NetCDF file `filename` with coordinate values `x_vals`, `y_vals` and `t_vals` or open an existing one.

    The format of the NetCDF file is chosen to exactly match the requirements of this program.
    """
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
                bdη_atts = Dict("units" => "m")

                x_dim = NcDim("x", x_vals, x_atts)
                y_dim = NcDim("y", y_vals, y_atts)
                t_dim = NcDim("time", t_vals, t_atts)

                b_var = NcVar("b", [x_dim, y_dim], atts=bdη_atts)

                u_var = NcVar("u", [x_dim, y_dim, t_dim], atts=uv_atts)
                v_var = NcVar("v", [x_dim, y_dim, t_dim], atts=uv_atts)
                d_var = NcVar("d", [x_dim, y_dim, t_dim], atts=bdη_atts)
                η_var = NcVar("eta", [x_dim, y_dim, t_dim], atts=bdη_atts)

                NetCDF.create(filename, b_var, u_var, v_var, d_var, η_var)
            else
                NetCDF.open(filename, mode=NC_WRITE)
            end

        NetCDF.nc_set_fill(nc.ncid, NetCDF.NC_NOFILL, Ptr{Int32}(0))
        return nc
    end
end
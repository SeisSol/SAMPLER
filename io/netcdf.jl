#=
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-07-13
=#

module NC
    using NetCDF
    using Main.Args

    export get_or_create_netcdf

    t_atts      = Dict("units" => "seconds")
    vel_atts    = Dict("units" => "m/s")
    dist_atts   = Dict("units" => "m")
    no_atts     = Dict()

    function get_netcdf(filename :: AbstractString) :: NcFile
        nc :: NcFile = NetCDF.open(filename, mode=NC_WRITE)
        NetCDF.nc_set_fill(nc.ncid, NetCDF.NC_NOFILL, Ptr{Int32}(0))
        return nc
    end

    """
    Create the NetCDF file `filename` with coordinate values `x_vals`, `y_vals` and `t_vals` or open an existing one.

    The format of the NetCDF file is chosen to exactly match the requirements of this program.
    """
    function create_netcdf(filename :: AbstractString, x_vals, y_vals, t_vals,
                           static_var_mappings  :: Array{Args.VarMapping, 1}, 
                           dynamic_var_mappings :: Array{Args.VarMapping, 1}) :: NcFile
        
        isfile(filename) && rm(filename)

        if isempty(x_vals) || isempty(y_vals) || isempty(t_vals) 
            throw(ArgumentError("Dimension values have to be passed when creating a NetCDF file!"))
        end

        function get_units(var_name:: String)
            if match(r"[uvw]", var_name) !== nothing
                return vel_atts
            elseif match(r"[UVW]", var_name) !== nothing
                return dist_atts
            else
                return no_atts
            end
        end

        x_dim = NcDim("x",    x_vals, dist_atts)
        y_dim = NcDim("y",    y_vals, dist_atts)
        t_dim = NcDim("time", t_vals, t_atts)

        static_vars  = []
        dynamic_vars = []

        for stat_mapping ∈ static_var_mappings
            append!(static_vars,  map(in_var -> NcVar(stat_mapping[in_var], [x_dim, y_dim],        atts=get_units(in_var)), collect(keys(stat_mapping))))
        end

        for dyn_mapping ∈ dynamic_var_mappings
            append!(dynamic_vars, map(in_var -> NcVar(dyn_mapping[in_var],  [x_dim, y_dim, t_dim], atts=get_units(in_var)), collect(keys(dyn_mapping))))
        end

        
        nc :: NcFile = NetCDF.create(filename, static_vars..., dynamic_vars...)

        NetCDF.nc_set_fill(nc.ncid, NetCDF.NC_NOFILL, Ptr{Int32}(0))
        return nc
    end
end

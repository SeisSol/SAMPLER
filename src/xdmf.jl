module XDMF
    using XMLDict
    using Mmap

    export grid_of, data_of

    struct XDMFFile
        xml         :: XMLDict.XMLDictElement
        base_path   :: AbstractString
        timesteps   :: Array{XMLDict.XMLDictElement,1}
    end

    function grid_of(file_path::AbstractString)
        xdmf = XDMFFile(file_path)

        # Memory-map geometry (points) and topology (tetrahedra) files.
        # Since the ids of the points referred to in topo_item start at 0 (and Julia counts from 1) we have to add 1.
        tetrahedra = mmap_data_item(xdmf.timesteps[1]["Topology"]["DataItem"], xdmf.base_path) .+ 1
        points     = mmap_data_item(xdmf.timesteps[1]["Geometry"]["DataItem"], xdmf.base_path)
        times = map(data_item -> parse(Float64, data_item["Time"][:Value]), xdmf.timesteps)

        return (tetrahedra, points, times)
    end

    function data_of(file_path::AbstractString, timestep::Int, vars...)
        xdmf = XDMFFile(file_path)
        timestep_attrs = xdmf.timesteps[timestep]["Attribute"]

        ret_array = Array{AbstractArray}(undef, length(vars))
        for (varid, var) ∈ enumerate(vars)
            for i ∈ 1:length(timestep_attrs)
                if timestep_attrs[i][:Name] == var
                    ret_array[varid] = mmap_hyperslab(timestep_attrs[i]["DataItem"], xdmf.base_path)
                    break
                end
            end
        end

        return Tuple(ret_array)
    end

    function XDMFFile(file_path::AbstractString)
        if !endswith(file_path, ".xmf") && !endswith(file_path, ".xdmf")
            throw(ArgumentError("The given file path has no ending associated with XDMF files (must be one of .xmf, .xdmf)"))
        end

        xml_file = open(file_path) do file; read(file, String) end
        xml_file = parse_xml(xml_file)

        # The root node "Xdmf" is assumed
        timesteps = xml_file["Domain"]["Grid"]["Grid"]
        num_timesteps = length(timesteps)

        if num_timesteps == 0
            error("No timesteps found in input file.")
        end

        xdmf_base_path = splitdir(file_path)[1]

        return XDMFFile(xml_file, xdmf_base_path, timesteps)
    end
        
    function mmap_hyperslab(data_item::XMLDict.XMLDictElement, base_path::AbstractString)
        #=
        <DataItem ItemType="HyperSlab" Dimensions="4561037">
            <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 0 1 1 1 4561037</DataItem>
            <DataItem NumberType="Float" Precision="8" Format="Binary" Dimensions="1 5242880">out_cell/mesh0/trunc_u.bin</DataItem>
        </DataItem>
        =#
        @assert data_item[:ItemType] == "HyperSlab"
        range_item = data_item["DataItem"][1]
        binary_item = data_item["DataItem"][2]

        @assert range_item[:Format] == "XML"
        @assert range_item[:Dimensions] == "3 2"
        @assert binary_item[:Format] == "Binary"

        target_dims = parse(UInt, data_item[:Dimensions])
        hyperslab_range = range_item[""]
        hyperslab_range = split(hyperslab_range, ' ')
        hyperslab_range = parse.(UInt, hyperslab_range)
        hyperslab_range = reshape(hyperslab_range, (2, 3))'
        @assert all(hyperslab_range[2,:] .== 1)

        hyperslab_range = [hyperslab_range[1,:]'; hyperslab_range[3,:]'] # Cut out stride row
        hyperslab_range[1,:] .+= 1 # Julia indices
        number_type = get_number_type(binary_item)

        file_range = binary_item[:Dimensions]
        file_range = split(file_range, ' ')
        file_range = parse.(UInt, file_range)
        file_range = tuple(reverse(file_range)...) # Column-major order!

        @assert file_range[2] == 1

        filename = binary_item[""]
        filename = joinpath(base_path, filename)
        whole_file = Mmap.mmap(filename, Array{number_type, 2}, file_range)
        return @view whole_file[hyperslab_range[1,2]:hyperslab_range[1,2]+hyperslab_range[2,2]-1,hyperslab_range[1,1]]
    end

    function mmap_data_item(data_item::XMLDict.XMLDictElement, base_path::AbstractString)
        number_type = get_number_type(data_item)        

        dimensions = data_item[:Dimensions]
        dimensions = split(dimensions, ' ')
        dimensions = parse.(UInt, dimensions)
        dimensions = tuple(reverse(dimensions)...) # Column-major order!

        filename = data_item[""]
        filename = joinpath(base_path, filename)
        return Mmap.mmap(filename, Array{number_type,length(dimensions)}, dimensions)
    end

    function get_number_type(data_item::XMLDict.XMLDictElement) :: Type
        bytes_of_precision = parse(UInt, data_item[:Precision])
        number_type_name = data_item[:NumberType]
        if bytes_of_precision == 8
            if number_type_name == "Int"
                return Int64
            elseif number_type_name == "UInt"
                return UInt64
            elseif number_type_name == "Float"
                return Float64
            else
                error("""Type "$number_type_name" of size $bytes_of_precision Bytes is not recognized by the parser.""")
            end
        elseif bytes_of_precision == 4
            if number_type_name == "Int"
                return Int32
            elseif number_type_name == "UInt"
                return UInt32
            elseif number_type_name == "Float"
                return Float32
            else
                error("""Type "$number_type_name" of size $bytes_of_precision Bytes is not recognized by the parser.""")
            end
        else
            error("""Type "$number_type_name" of size $bytes_of_precision Bytes is not recognized by the parser.""")
        end
    end
end
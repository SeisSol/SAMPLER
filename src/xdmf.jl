module XDMF
    using XMLDict
    using Mmap
    using Base.Threads
    using LinearAlgebra
    using Test
    include("quadtrees.jl")

    export XDMFFile, PolygonBounds, mmap_data_item

    struct SparseSamplingRaster
        raster_cells::Dict{Tuple{Int32, Int32}, Array{Int32,1}}
        x_min::Float64; x_max::Float64; dx::Float64
        y_min::Float64; y_max::Float64; dy::Float64
    end

    struct DenseSamplingRaster
        raster_cells::Array{Array{Int32,1}, 2}
        x_min::Float64; x_max::Float64; dx::Float64
        y_min::Float64; y_max::Float64; dy::Float64
    end

    struct PolygonBounds
        x_min::Float64; x_max::Float64
        y_min::Float64; y_max::Float64
    end

    mutable struct T0T1
        t0 :: Float64
        t1 :: Float64
    end

    function XDMFFile(file_path::AbstractString, num_sampling_columns::Tuple{Integer,Integer})
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

        # Memory-map geometry (points) and topology (tetrahedra) files.
        # Since the ids of the points referred to in topo_item start at 0 (and Julia counts from 1) we have to add 1.
        tetrahedra = mmap_data_item(timesteps[1]["Topology"]["DataItem"], xdmf_base_path) .+ 1
        points     = mmap_data_item(timesteps[1]["Geometry"]["DataItem"], xdmf_base_path)

        domain_x = (minimum(points[1,:]), maximum(points[1,:]))
        domain_y = (minimum(points[2,:]), maximum(points[2,:]))

        println("Domain is $domain_x×$domain_y.")

        sampling_dx = (domain_x[end] - domain_x[1]) / num_sampling_columns[1]
        sampling_dy = (domain_y[end] - domain_y[1]) / num_sampling_columns[2]

        # shape(tetrahedra) = 4×n; iterate over n
        region = Quadtrees.Quad(domain_x[1], domain_y[1], domain_x[2] - domain_x[1], domain_y[2] - domain_y[1])
        qt = Quadtrees.Quadtree(region, tetrahedra, points)

        # return XDMFFile()
    end

    function is_point_in_tetrahedron(p::AbstractArray{Float64, 1}, tetrahedron::AbstractArray{Float64,2}) :: Bool
        for i ∈ 1:4
            # Exclude point i from the tetrahedron and consider the plane defined by the other 3 points
            # If and only if (x,y,z) lies on the "positive" side (in normal direction) of all 4 planes, 
            # it is inside the tetrahedron
            p_excl = tetrahedron[i,:]
            p1 = tetrahedron[(i%4)+1,:]
            p2 = tetrahedron[((i+1)%4)+1,:]
            p3 = tetrahedron[((i+2)%4)+1,:]

            # Calculate Hesse normal form of the plane defined by p1, p2, p3
            # Normal vector
            n = normalize(cross(p2-p1, p3-p1))

            # n1x1 + n2x2 + n3x3 + d = 0; solve for d
            d = -dot(p1, n)

            # If p_excl is ALWAYS on the "positive" side of the plane. 
            # So, if p_excl and p lie on the same side, p is inside the tetrahedron
            dist_p = dot(p, n) + d
            dist_p_excl = dot(p_excl, n) + d

            # If signs of distance vectors do not match, p and p_excl are on different sides of the plane (reject)
            if ((dist_p <= 0 && dist_p_excl > 0) || (dist_p > 0 && dist_p_excl <= 0))
                return false
            end
        end

        return true
    end

    function mmap_data_item(data_item::XMLDict.XMLDictElement, base_path::AbstractString)
        bytes_of_precision = parse(UInt, data_item[:Precision])
        number_type_name = data_item[:NumberType]
        number_type =
            if bytes_of_precision == 8
                if number_type_name == "Int"
                    Int64
                elseif number_type_name == "UInt"
                    UInt64
                elseif number_type_name == "Float"
                    Float64
                else
                    error("""Type "$number_type_name" of size $bytes_of_precision Bytes is not recognized by the parser.""")
                end
            elseif bytes_of_precision == 4
                if number_type_name == "Int"
                    Int32
                elseif number_type_name == "UInt"
                    UInt32
                elseif number_type_name == "Float"
                    Float32
                else
                    error("""Type "$number_type_name" of size $bytes_of_precision Bytes is not recognized by the parser.""")
                end
            else
                error("""Type "$number_type_name" of size $bytes_of_precision Bytes is not recognized by the parser.""")
            end

        dimensions = data_item[:Dimensions]
        dimensions = split(dimensions, ' ')
        dimensions = parse.(UInt, dimensions)
        dimensions = tuple(reverse(dimensions)...) # Why does Julia fill matrices column-wise!? 😭

        filename = data_item[""] # No joke, that's how the XMLDict library stores the text in an XML tag 😂
        filename = joinpath(base_path, filename)
        return Mmap.mmap(filename, Array{number_type,length(dimensions)}, dimensions)
    end
end
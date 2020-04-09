module XDMF
    using XMLDict
    using Mmap
    using Base.Threads
    using LinearAlgebra
    using Test
    using Traceur
    include("quadtrees.jl")

    export XDMFFile

    struct SamplingRaster
        raster_cells::Array{Array{NTuple{2, Int32}, 1},2}
        x_min::Float64; x_max::Float64; dx::Float64
        y_min::Float64; y_max::Float64; dy::Float64
    end

    mutable struct ThreadData
        cells_processed     :: Int32
        p                   :: Array{Float64,1}
        tet_dict            :: Dict{Int32, Int32}
        tetrahedra          :: Array{Quadtrees.Tetrahedron, 1}
    end

    function XDMFFile(file_path::AbstractString, num_sampling_columns::NTuple{2,Integer})
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
        domain_z = (minimum(points[3,:]), maximum(points[3,:]))

        println("Domain is $domain_x × $domain_y × $domain_z.")

        sampling_dx = (domain_x[end] - domain_x[1]) / num_sampling_columns[1]
        sampling_dy = (domain_y[end] - domain_y[1]) / num_sampling_columns[2]
        sampling_dz = (domain_z[end] - domain_z[1]) / 40

        # shape(tetrahedra) = 4×n; iterate over n
        region = Quadtrees.Quad(domain_x[1], domain_y[1], domain_x[2] - domain_x[1], domain_y[2] - domain_y[1])
        qt = Quadtrees.Quadtree(region, tetrahedra, points)

        tetrahedra = nothing
        points = nothing

        println("Allocating cell_array...")

        cell_array = [Array{NTuple{2, Int32}, 1}() for i ∈ 1:num_sampling_columns[2], j ∈ 1:num_sampling_columns[1]]
        raster = SamplingRaster(cell_array, 
                                domain_x[1], domain_x[2], sampling_dx,
                                domain_y[1], domain_y[2], sampling_dy)

        println("Done allocating cell_array.")

        print_lock = SpinLock()
        sc_n = num_sampling_columns[1] * num_sampling_columns[2] / nthreads()

        thread_contexts = [ThreadData(0, [0,0,0], Dict{Int32, Int32}(), Array{Quadtrees.Tetrahedron, 1}(undef, 2048)) for id ∈ 1:nthreads()]

        Threads.@threads for x ∈ 1:num_sampling_columns[1]
            ctxt = thread_contexts[threadid()]

            ctxt.p[1] = domain_x[1] + sampling_dx * (x - 1)

            for y ∈ 1:num_sampling_columns[2]
                if ctxt.cells_processed % 1024 == 0
                    lock(print_lock)
                    Core.println("Processing cell (y=",y," x=",x,")... ",ctxt.cells_processed/sc_n*100,"% done.") 
                    unlock(print_lock)
                end

                ctxt.p[2] = domain_y[1] + sampling_dy * (y - 1)

                Quadtrees.get_tetrahedra_at!(qt, ctxt.p[1], ctxt.p[2], ctxt.tetrahedra)
                cell :: AbstractArray{NTuple{2, Int32}, 1} = cell_array[y, x]
                empty!(ctxt.tet_dict)

                for z ∈ 1:40
                    ctxt.p[3] = domain_z[1] + sampling_dz * (z - 1)
                    for tet ∈ ctxt.tetrahedra
                        if is_point_in_tetrahedron(ctxt.p, tet.points)
                            if !haskey(ctxt.tet_dict, tet.id); ctxt.tet_dict[tet.id] = 1
                            else ctxt.tet_dict[tet.id] += 1 end
                        end
                    end
                end

                append!(cell, [(pair.first, pair.second) for pair :: Pair{Int32, Int32} ∈ pairs(ctxt.tet_dict)])
                ctxt.cells_processed += 1
            end
        end

        show(IOContext(stdout, :limit=>true), "text/plain", cell_array)
        # return XDMFFile()
    end

    function is_point_in_tetrahedron(p::AbstractArray{Float64, 1}, tetrahedron::AbstractArray{Float64,2}) :: Bool
        @views points = (tetrahedron[1,:], tetrahedron[2,:], tetrahedron[3,:], tetrahedron[4,:])
        n = [0., 0., 0.]

        function cross!(a, b, ret)
            ret[1] = (a[2]*b[3] - a[3]*b[2])
            ret[2] = (a[3]*b[1] - a[1]*b[3])
            ret[3] = (a[1]*b[2] - a[2]*b[1])
        end
        
        for i ∈ 1:4
            # Exclude point i from the tetrahedron and consider the plane defined by the other 3 points
            # If and only if (x,y,z) lies on the "positive" side (in normal direction) of all 4 planes, 
            # it is inside the tetrahedron
            p_excl = points[i]
            p1 = points[(i%4)+1]
            p2 = points[((i+1)%4)+1]
            p3 = points[((i+2)%4)+1]

            # Calculate Hesse normal form of the plane defined by p1, p2, p3
            # Normal vector
            cross!(p2-p1, p3-p1, n)
            normalize!(n)

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
module XDMF
    using XMLDict
    using Mmap
    using Base.Threads
    using Traceur

    export XDMFFile, PolygonBounds, rasterize_convex_hull!, mmap_data_item, calculate_convex_hull_points, calculate_overlap!

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

        rasters = [
            SparseSamplingRaster(Dict{Tuple{Int32, Int32}, Array{Int32, 1}}(), 
                           domain_x[1], domain_x[end], sampling_dx,
                           domain_y[1], domain_y[end], sampling_dy)
            for i ∈ 1:nthreads()
        ]

        total_tetrahedra_processed = Threads.Atomic{Int}(0)
        last_printed_processed = 0
        
        # shape(tetrahedra) = 4×n; iterate over n
        num_tetrahedra = size(tetrahedra, 2)

        Threads.@threads for i ∈ 1:num_tetrahedra
            if threadid() == 2
                processed = total_tetrahedra_processed[]
                if processed - last_printed_processed ≥ 10000
                    println("Processed $processed of $num_tetrahedra tetrahedra... $(processed/num_tetrahedra*100)% done.")
                    flush(stdout)
                    last_printed_processed = processed
                end
            end

            tetrahedron = tetrahedra[:, i] # shape(tetrahedron) = 4×1
            tetrahedron = [points[dim, point_id] for point_id in tetrahedron, dim in 1:3] # shape(tetrahedron) = 4×3

            calculate_overlap!(tetrahedron, rasters[threadid()], i)
            atomic_add!(total_tetrahedra_processed, 1)
        end

        # Collect sampling data from threads
        raster = DenseSamplingRaster([Array{Int32,1}() for i ∈ 1:num_sampling_columns[1], j ∈ 1:num_sampling_columns[2]],
                                     domain_x[1], domain_x[end], sampling_dx,
                                     domain_y[1], domain_y[end], sampling_dy)

        for i ∈ 1:length(rasters)
            r = rasters[i]
            for ((y, x), ls) ∈ pairs(r.raster_cells)
                append!(raster.raster_cells[y, x], ls)
            end

            rasters[i] = nothing
        end

        show(IOContext(stdout, :limit=>true), "text/plain", raster.raster_cells)
        # return XDMFFile()
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

    function calculate_overlap!(tetrahedron::AbstractArray{Float64,2}, raster::SparseSamplingRaster, tetrahedron_id::Int)
        px_min = minimum(tetrahedron[:,1]); px_max = maximum(tetrahedron[:,1])
        py_min = minimum(tetrahedron[:,2]); py_max = maximum(tetrahedron[:,2])

        # A list of either 1 or 2 disjoint triangles that make up the tetrahedron's complex hull
        # when projection it onto the XY-Plane
        convex_hull = calculate_convex_hull_points(tetrahedron)
        rasterize_convex_hull!(convex_hull, PolygonBounds(px_min, px_max, py_min, py_max), raster, tetrahedron_id)
    end

    function calculate_convex_hull_points(tetrahedron::AbstractArray{Float64,2}) :: AbstractArray{Float64,2}
        # Calculate convex hull (in top-down view, 2D) with Graham's Algorithm
        # The anchor point p0 is the point with the lowest y-coordinate.
        # Should there be >1 lowest points, the one with the lowest x-coordinate of them is chosen.
        # Therefore, sort by x-coord first and then take the first point with the lowest y-coordinate
        tetrahedron = sortslices(tetrahedron, dims = 1, lt = (a, b)->a[1] < b[1])
        index_y_min = argmin(tetrahedron[:,2])
        p0 = tetrahedron[index_y_min,:]

        function p1_clockwise_of_p2(p1::AbstractArray{Float64,1}, p2::AbstractArray{Float64,1}) :: Bool
            # Rotate p1 90° CCW around origin: p1' := (-p1y, p1x)
            # The scalar product p1'*p2 is then 0 if they are aligned, > 0 if p1 is CW from p2 and < 0 otherwise
            -p1[2]*p2[1] + p1[1]*p2[2] > 0
        end

        remaining = [tetrahedron[1:index_y_min - 1,:]; tetrahedron[index_y_min + 1:end,:]] # cut row of p0 out of list
        # sort by polar angle relative to x-axis parallel through p0
        remaining = sortslices(remaining, dims = 1, lt = p1_clockwise_of_p2)

        # if the angle curls clockwise at any point (turns out, when 3 points sorted by angle are left, 
        # that can only be the angle below), remove that point from the convex hull
        if p1_clockwise_of_p2(remaining[1,:] .- remaining[2,:], remaining[3,:] .- remaining[2,:])
            remaining = [transpose(remaining[1,:]); transpose(remaining[3,:])]
        end

        # This is our convex hull in counter-clockwise order!
        return [transpose(p0); remaining]
    end

    function rasterize_convex_hull!(convex_hull::AbstractArray{Float64,2}, bounds::PolygonBounds, raster::SparseSamplingRaster, tetrahedron_id::Int)
        # This is the bounding box of all grid columns that could be affected by the tetrahedron
        px_min = bounds.x_min; px_max = bounds.x_max
        py_min = bounds.y_min; py_max = bounds.y_max

        # Index of the row/column that contains the respective coordinate 'num', in 'raster'
        idx(num::Float64, delta::Float64, offset::Float64)::Int = floor(Int, (num - offset) / delta) + 1

        # The index offset (relative to 'raster') the 'overlap_grid' below starts at.
        overlap_start_idx_x = idx(px_min, raster.dx, raster.x_min)
        overlap_start_idx_y = idx(py_min, raster.dy, raster.y_min)
        overlap_end_idx_x = idx(px_max, raster.dx, raster.x_min)
        overlap_end_idx_y = idx(py_max, raster.dy, raster.y_min)

        nx = overlap_end_idx_x - overlap_start_idx_x
        ny = overlap_end_idx_y - overlap_start_idx_y

        overlap_grid = zeros(Bool, (ny, nx))

        # This function returns the index in the above grid obtained from float coordinates lying inside it
        # +1 as Julia is 1-indexed
        overlap_idx_x(x::Float64) = idx(x, raster.dx, raster.x_min) - overlap_start_idx_x + 1
        overlap_idx_y(y::Float64) = idx(y, raster.dy, raster.y_min) - overlap_start_idx_y + 1

        overlap_pos_x(idx_x::Number) = idx_x * raster.dx + raster.x_min
        overlap_pos_y(idx_y::Number) = idx_y * raster.dy + raster.y_min

        # For each cell (x,y) in the domain of overlap_grid, check overlap with polygon
        Threads.@threads for idx_x ∈ overlap_start_idx_x:overlap_end_idx_x-1
            left = overlap_pos_x(idx_x); right = overlap_pos_x(idx_x + 1)
        
            for idx_y ∈ overlap_start_idx_y:overlap_end_idx_y-1
                bottom = overlap_pos_y(idx_y); top = overlap_pos_y(idx_y + 1)

                # For each line between two individual points
                n_points = size(convex_hull, 1)
                for i ∈ 1:n_points
                    (x0, y0) = Tuple(convex_hull[i,1:2]) # Only x- and y-coordinate are of importance
                    (x1, y1) = Tuple(convex_hull[i % n_points + 1,1:2])

                    # Mark cell according to Liang-Barsky clipping algorithm
                    # Adapted from the Algorithm on page 7 of their "A New Concept and Method for Line Clipping" paper from 1984
                    function clipt!(p::Float64, q::Float64, t0t1::T0T1)
                        if p < 0
                            r = q/p
                            if r > t0t1.t1; return false
                            elseif r > t0t1.t0; t0t1.t0 = r end
                        elseif p > 0
                            r = q/p
                            if r < t0t1.t0; return false
                            elseif r < t0t1.t1; t0t1.t1 = r end
                        else # p = 0
                            if q < 0; return false end
                        end

                        return true
                    end

                    function line_intersects_square()
                        t0t1 = T0T1(0., 1.)
                        Δx = x1 - x0
                        if clipt!(-Δx, x0 - left, t0t1) && clipt!(Δx, right - x0, t0t1)
                            Δy = y1 - y0
                            return clipt!(-Δy, y0 - bottom, t0t1) && clipt!(Δy, top - y0, t0t1)
                        else
                            return false
                        end
                    end

                    if line_intersects_square()
                        overlap_grid[idx_y-overlap_start_idx_y+1, idx_x-overlap_start_idx_x+1] = true
                    end
                end
            end
        end

        # Flood-fill the convex hull, the inside of which is now fully enclosed by the marked cells on its borders
        for x ∈ 1:size(overlap_grid, 2)
            y0 = 1; y1 = size(overlap_grid, 1)

            # Find starting and ending marked cell in the current column. If found, fill all between
            while y0 < y1 && !overlap_grid[y0, x]; y0 += 1 end
            while y0 < y1 && !overlap_grid[y1, x]; y1 -= 1 end
            # Start on the cell AFTER the starting marked cell (no need to re-mark it)
            while (y0 += 1) < y1; overlap_grid[y0, x] = true end # Fill cells BETWEEN (excl. borders) starting and ending marked cell (if applicable) 
        end

        for x ∈ 1:size(overlap_grid, 2), y ∈ 1:size(overlap_grid, 1)
            if overlap_grid[y, x]
                key = (overlap_start_idx_y+y-1, overlap_start_idx_x+x-1)
                if !haskey(raster.raster_cells, key)
                    raster.raster_cells[key] = [tetrahedron_id]
                else
                    push!(raster.raster_cells[key], tetrahedron_id)
                end
            end
        end
    end
end
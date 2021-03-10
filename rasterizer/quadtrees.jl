#=
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-04-01
=#

module Quadtrees
    export Quadtree, quadtree_insert!, Quad, Tetrahedron, get_tetrahedra_at!

    struct Quad
        x        :: Float64
        y        :: Float64
        width    :: Float64
        height   :: Float64
    end

    struct Tetrahedron
        id      :: Int32
        aabb    :: NTuple{4, Float64}
        points  :: AbstractArray{Float64, 2}
    end

    mutable struct Quadtree
        region   :: Quad
        contents :: Array{Tetrahedron,1}
        xn_yn    :: Union{Quadtree, Nothing}
        xp_yn    :: Union{Quadtree, Nothing}
        xn_yp    :: Union{Quadtree, Nothing}
        xp_yp    :: Union{Quadtree, Nothing}
    end

    @enum Quadrant center xn_yn xp_yn xn_yp xp_yp

    function Quadtree(region::Quad, tetrahedra::AbstractArray{Int64, 2}, points::AbstractArray{Float64, 2})
        qt = Quadtree(deepcopy(region), 
                      Array{Tetrahedron,1}(), 
                      nothing, nothing, nothing, nothing)

        num_tetrahedra = size(tetrahedra, 2)
        for tet_id :: Int32 ∈ 1:num_tetrahedra
            if tet_id % 10000 == 0; println("Progress: tet $tet_id of $num_tetrahedra, $(tet_id/num_tetrahedra*100)% done") end

            tetrahedron = tetrahedra[:, tet_id]
            tet_points = [points[dim, point_id] for point_id in tetrahedron, dim in 1:3] # shape(t_points) = (4×3)
            tet_aabb = (minimum(tet_points[:, 1]), minimum(tet_points[:, 2]), 
                        maximum(tet_points[:, 1]), maximum(tet_points[:, 2]))
            quadtree_insert!(qt, Tetrahedron(tet_id, tet_aabb, tet_points))
        end

        return qt
    end

    function quadtree_insert!(qt::Quadtree, tet::Tetrahedron)
        q = get_quadrant(qt.region, tet.aabb)
        if q == center 
            push!(qt.contents, tet)
            return
        end

        current_elem = qt[q]

        if isnothing(current_elem)
            new_width = qt.region.width / 2
            new_height = qt.region.height / 2

            # If the x-component of our quadrant is negative, copy the left border of the current
            # tree's region. Otherwise, add half of the current tree's width to it. Same for y.
            new_x = (q == xn_yn || q == xn_yp) ? qt.region.x : qt.region.x + new_width
            new_y = (q == xn_yn || q == xp_yn) ? qt.region.y : qt.region.y + new_height

            new_region = Quad(new_x, new_y, new_width, new_height)

            new_tree = Quadtree(new_region, 
                                Array{Tetrahedron,1}(), 
                                nothing, nothing, nothing, nothing)

            qt[q] = new_tree
            quadtree_insert!(new_tree, tet)
            quadtree_insert!(new_tree, tet)
        else # current_elem isa Quadtree
            quadtree_insert!(current_elem, tet)
        end
    end

    function get_tetrahedra_at!(qt::Quadtree, x::Float64, y::Float64, ret_tetrahedra::Array{Tetrahedron, 1})
        empty!(ret_tetrahedra)
        current_node :: Union{Quadtree, Nothing} = qt
        while !isnothing(current_node)
            for tet :: Tetrahedron ∈ current_node.contents
                if x ≥ tet.aabb[1] && x ≤ tet.aabb[3] && y ≥ tet.aabb[2] && y ≤ tet.aabb[4]
                    push!(ret_tetrahedra, tet)
                end
            end

            quadrant :: Quadrant = get_quadrant(current_node.region, x, y)
            current_node = current_node[quadrant]
        end
    end

    function get_quadrant(region::Quad, bounds::NTuple{4, Float64}) :: Quadrant
        # For numerical accuracy, coordinates are subtracted before comparing position within region
        # (coordinates have a larger domain than coordinates within a region)

        # x_max is on the left of the center?
        is_left = bounds[3] - region.x <= region.width * .5
        # y_max is at the bottom of the center?
        is_btm = bounds[4] - region.y <= region.height * .5

        if !is_left
            is_right = bounds[1] - region.x > region.width * .5
            if !is_right return center end
        end

        if !is_btm
            is_top = bounds[2] - region.y > region.height * .5
            if !is_top return center end
        end

        if is_left  && is_btm  return xn_yn end
        if is_left  && !is_btm return xn_yp end
        if !is_left && is_btm  return xp_yn end
        if !is_left && !is_btm return xp_yp end

        return xn_yn # NEVER happens but compiler is now happy.
    end

    function get_quadrant(region::Quad, x::Float64, y::Float64) :: Quadrant
        is_left = x - region.x <= region.width * .5
        is_btm = y - region.y <= region.height * .5

        if is_left  && is_btm  return xn_yn end
        if is_left  && !is_btm return xn_yp end
        if !is_left && is_btm  return xp_yn end
        if !is_left && !is_btm return xp_yp end

        return xn_yn # NEVER happens but compiler is now happy.
    end

    function Base.getindex(collection::Quadtree, key::Quadrant) :: Union{Quadtree, Nothing}
        if key == xn_yn return collection.xn_yn end
        if key == xn_yp return collection.xn_yp end
        if key == xp_yn return collection.xp_yn end
        if key == xp_yp return collection.xp_yp end

        throw(ArgumentError("Illegal key: $key"))
    end

    function Base.setindex!(collection::Quadtree, value::Union{Quadtree, Nothing}, key::Quadrant)
        if key == xn_yn collection.xn_yn = value
        elseif key == xn_yp collection.xn_yp = value
        elseif key == xp_yn collection.xp_yn = value
        elseif key == xp_yp collection.xp_yp = value
        else throw(ArgumentError("Illegal key: $key")) end
    end
end

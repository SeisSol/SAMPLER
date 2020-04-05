#=
- Julia version: 1.4.0
- Author: Maximilian Schmeller
- Date: 2020-04-01
=#

module Quadtrees
    export Point, Quadtree, quadtree_insert!, Quad

    struct Point{T}
        x       :: Float64
        y       :: Float64
        content :: T
    end

    struct Quad
        x       :: Float64
        y       :: Float64
        width   :: Float64
        height  :: Float64
    end

    mutable struct Quadtree{T}
        region  :: Quad
        xn_yn   :: Union{Quadtree{T}, Point{T}, Nothing}
        xp_yn   :: Union{Quadtree{T}, Point{T}, Nothing}
        xn_yp   :: Union{Quadtree{T}, Point{T}, Nothing}
        xp_yp   :: Union{Quadtree{T}, Point{T}, Nothing}
    end

    @enum Quadrant xn_yn xp_yn xn_yp xp_yp

    function Quadtree(region::Quad, points::AbstractArray{Point{T}}) where T
        qt = Quadtree{T}(deepcopy(region), nothing, nothing, nothing, nothing)

        for point in points
            quadtree_insert!(qt, point)
        end
    end

    function quadtree_insert!(qt::Quadtree{T}, point::Point{T}) where T
        q = get_quadrant(qt.region, point.x, point.y)
        current_elem = qt[q]

        if isnothing(current_elem)
            qt[q] = point
        elseif current_elem isa Quadtree
            quadtree_insert!(current_elem, point)
        else # current_elem isa Point
            new_width = qt.region.width / 2
            new_height = qt.region.height / 2

            # If the x-component of our quadrant is negative, copy the left border of the current
            # tree's region. Otherwise, add half of the current tree's width to it. Same for y.
            new_x = (q == xn_yn || q == xn_yp) ? qt.region.x : qt.region.x + new_width
            new_y = (q == xn_yn || q == xp_yn) ? qt.region.y : qt.region.y + new_height

            new_region = Quad(new_x, new_y, new_width, new_height)

            new_tree = Quadtree{T}(new_region, nothing, nothing, nothing, nothing)
            qt[q] = new_tree
            quadtree_insert!(new_tree, current_elem)
            quadtree_insert!(new_tree, point)
        end
    end

    function get_quadrant(region::Quad, x::Real, y::Real)
        # For numerical accuracy, coordinates are subtracted before comparing position within region
        # (coordinates have a larger domain than coordinates within a region)
        is_left = x - region.x <= region.width / 2
        is_btm = y - region.y <= region.height / 2

        if is_left  && is_btm  return xn_yn end
        if is_left  && !is_btm return xn_yp end
        if !is_left && is_btm  return xp_yn end
        if !is_left && !is_btm return xp_yp end
    end

    function Base.getindex(collection::Quadtree, key::Quadrant)
        if key == xn_yn return collection.xn_yn end
        if key == xn_yp return collection.xn_yp end
        if key == xp_yn return collection.xp_yn end
        if key == xp_yp return collection.xp_yp end
    end

    function Base.setindex!(collection::Quadtree{T}, value::Union{Quadtree{T}, Point{T}, Nothing}, key::Quadrant) where T
        if key == xn_yn collection.xn_yn = value end
        if key == xn_yp collection.xn_yp = value end
        if key == xp_yn collection.xp_yn = value end
        if key == xp_yp collection.xp_yp = value end
    end
end

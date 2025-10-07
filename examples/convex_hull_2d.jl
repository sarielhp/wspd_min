## Silly julia code to compute convex hull in 2d

@inline function  is_leq_turn( p::Point2F, q::Point2F, r::Point2F )
    q_x = q[1] - p[1];
    q_y = q[2] - p[2];
    r_x = r[1] - p[1];
    r_y = r[2] - p[2];

    return   ( q_x * r_y - q_y * r_x ) >=  0.0;
end

function convex_hull_2d( _points::Vector{Point2F} )
    n = length( _points )
    if n <= 3
        return unique(_points) # Handle small cases and duplicates
    end

    # Step 1: Sort points by x-coordinate, with y as a tie-breaker
    points = sort(_points, by = p -> (p[1], p[2]))
    unique!( points ) # remove duplicates
    n = length( points )
    
    # Step 2: Build the lower hull
    lower_hull = Point[]
    for p in points
        while length(lower_hull) >= 2  &&  is_leq_turn(lower_hull[end-1], lower_hull[end], p)
            pop!(lower_hull)
        end
        push!(lower_hull, p)
    end

    # Step 3: Build the upper hull
    upper_hull = Point[]
    for i in n:-1:1
        p = points[i]
        while length(upper_hull) >= 2 && is_leq_turn(upper_hull[end-1], upper_hull[end], p)
            pop!(upper_hull)
        end
        push!(upper_hull, p)
    end

    # Step 4: Combine the two hulls, removing duplicates at the start and end
    # The last point of the lower hull and the last point of the upper hull are the same
    # as the first point of the upper and lower hulls, respectively.
    pop!(lower_hull)
    println( lower_hull )
    pop!(upper_hull)
    
    return [lower_hull; upper_hull]
end

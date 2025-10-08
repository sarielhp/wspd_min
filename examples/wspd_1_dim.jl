#! /bin/env  julial

push!(LOAD_PATH, pwd()*"/src/")
#push!(LOAD_PATH, pwd()*"/src/cg/")

using FrechetDist
using FrechetDist.cg
using FrechetDist.cg.polygon
using FrechetDist.cg.bbox
using FrechetDist.cg.point
using VirtArray
using Printf
using DataFrames
using PrettyTables
using Random

using JuMP, Gurobi
#using HiGHS, GLPK, JuMP, Gurobi
include( "graphics.jl" )

import BBT
import WSPD
using Cairo, Colors

include( "convex_hull_2d.jl" )

#########################################################################################

Point1F = Point{1,Float64}
Polygon1F = Polygon{1,Float64}

################################################################################
# stuff to add bbox.jl
################################################################################


function  BBox{D,T}( p::Point{D,T}, q::Point{D,T} ) where{D,T}
    bb = BBox{D,T}()
    init( bb, p, q );
    return  bb
end

function   Base.:*(x_r::Tuple{T,T}, y_r::Tuple{T, T}) where{T}
    return  BBox{2,T}( Point{2,T}( x_r[ 1 ], y_r[1] ),
                       Point{2,T}( x_r[ 2 ], y_r[2] ) )
end

function   Base.:*(b::BBox{D,T}, rng::Tuple{T, T})::BBox{D+1,T} where{D,T}
    mini = Point{D+1,T}( bottom_left( b )..., rng[1] )
    maxi = Point{D+1,T}( top_right( b )..., rng[2] )
    return  BBox{D+1,T}( mini, maxi )
end

function  vertices( bb::BBox2F )
    return  [ bottom_left( bb ), bottom_right( bb ),
              top_left( bb ), top_right( bb ) ]
end


function  Base.isempty( b::BBox{D,T} )::Bool  where{D,T}
    if  ! b.f_init
        return  true 
    end
    for  i ∈ 1:D
        if  d_min( b, i ) > d_max( b, i )
            return  true
        end
    end
    return  false
end


function  issubset( b::BBox{D,T}, c::BBox{D,T} )::Bool where {D,T}
    #println( "issubset" )
    if  isempty( b )
        return  true;
    end
    if   isempty( c )
        return  false;
    end
    for  i ∈ 1:D
        if  d_min( b, i ) < d_min( c, i )
            return  false
        end
        if  d_max( b, i ) > d_max( c, i )
            return  false
        end
        #@assert( d_min(c,i) <= d_min(b, i ) <= d_max(b, i) <= d_max( c, i ) )
    end

    return  true
end

function  ⊆( b::BBox{D,T}, c::BBox{D,T} )::Bool where {D,T}
    return  issubset( b, c )
end

function  get_dim( bb, d )
    return  ( bb.mini[ d ], bb.maxi[ d ] )
end

function   compatible( b::BBox{D,T}, c::BBox{D,T} )::Tuple{Bool,Int}  where {D,T}
    diff = 0
    dim = 0
    for i ∈ 1:D
        if  get_dim( b, i ) != get_dim(c, i )
            if  ( ( d_min(b,i) != d_max(c,i ) )
                  &&  ( d_max(b,i) != d_min(c,i ) ) )
                return  false, 0
            end
            diff += 1
            dim = i
            if  diff > 1
                return  false, 0
            end
        end
    end

    return  true, dim
end

function  merge( bb::BBox{D,T}, bbo::BBox{D,T} ) where {D,T}
    out = BBox{D,T}()
    for  b ∈ [bb, bbo ]
        bound( out, bottom_left( b ) )
        bound( out, top_right( b ) )
    end
    return   out
end

function  is_intersect( b::BBox{D,T}, c::BBox{D,T} )::Bool where {D,T}
    for i ∈ 1:D
        if  ( ( d_max( b, i ) <= d_min( c, i ) )
              ||  ( d_min( b, i ) >= d_max( c, i ) ) )
            return  false
        end
    end
    return  true;    
end

################################################################################
################################################################################

# The anchor stuff

function  energy( p::Point2F )::Float64
    @assert( p[2] >= p[ 1 ] )
    return  p[ 2 ] - p[ 1 ]
end

function   energy( bb::BBox2F )
    return  energy(  bottom_left( bb ) )
end


function  sq_from_anchor( p::Point2F, ε )
    ℓ = ε * energy( p )
    
    p_min = Point2F( p[ 1 ] - ℓ, p[2] )
    p_max = Point2F( p[ 1 ], p[2] + ℓ  )

    return  BBox2F_init( p_min, p_max )
end

function  sq_from_top_left( p::Point2F, ε )
    x, y = p[ 1 ], p[2] 
    new_e = energy( p )  / ( 1 + 2.0 * ε ) 
    ℓ = ( energy( p ) - new_e ) / 2.0

    new_anc = Point2F( x + ℓ, y - ℓ )
    return  sq_from_anchor( new_anc, ε )
end

function  sq_from_bottom_left( p::Point2F, ε )
    x, y = p[ 1 ], p[2] 
    new_e = energy( p )  / ( 1.0 + ε ) 
    ℓ = energy( p ) - new_e

    new_anc = Point2F( x + ℓ, y )
    return  sq_from_anchor( new_anc, ε )
end


function  sq_from_top_right( p::Point2F, ε )
    x, y = p[ 1 ], p[2] 
    new_e = energy( p )  / ( 1.0 + ε ) 
    ℓ = energy( p ) - new_e

    new_anc = Point2F( x, y - ℓ )
    return  sq_from_anchor( new_anc, ε )
end


function  cover_polygon( p::Point2F, ε, boxes::Vector{BBox2F} )
    top_right = sq_from_bottom_left( p, ε )
    top_left = sq_from_anchor( p, ε )
    bottom_right = sq_from_top_left( p, ε )
    bottom_left = sq_from_top_right( p, ε )

    P = Polygon2F()
    push!( P,  bbox.top_right( top_right ),
               bbox.top_right( top_left ),
               bbox.top_left( top_left ),
               bbox.bottom_left( top_left ),
               bbox.bottom_left( bottom_left ),
               bbox.bottom_right( bottom_left ),
               bbox.bottom_right( bottom_right ),
               bbox.bottom_right( top_right ) )
    push!( boxes, top_right, top_left, bottom_right, bottom_left )
    return  P
end


function  cover_for_point( p, ε )
    BS = Vector{BBox2F}()
    
    sq_anc   = sq_from_anchor( p, ε )
    sq_below = sq_from_top_right( p, ε )

    q = bottom_left( sq_below )
    sq_b_l = sq_from_anchor( q, ε )
    
    push!( BS, sq_anc, sq_below, sq_b_l )

#    exit( -1 )
    return  BS
end


##########################################################################################
# the set-system stuff was written by Gemini...
#
mutable struct SetSystem
    n_items::Int
    m_sets::Int
    sets::Matrix{Bool}

    # Constructor
    function SetSystem(n_items::Int, m_sets::Int)
        # Check for valid input
        if n_items <= 0 || m_sets <= 0
            error("Number of items and sets must be positive integers.")
        end
        # Initialize the boolean matrix with all false values
        sets = fill(false, n_items, m_sets)
        new(n_items, m_sets, sets)
    end
end

"""
    Base.setindex!(ss::SetSystem, value::Bool, set_index::Int, item_index::Int)

Sets whether a specific item is in a specific set within the SetSystem
using the `ss[set_index, item_index] = value` syntax.

# Arguments
- `ss::SetSystem`: The SetSystem instance to modify.
- `value::Bool`: A boolean value. If `true`, the item is added. If `false`, it is removed.
- `set_index::Int`: The index of the set to modify (1-based).
- `item_index::Int`: The index of the item to add or remove (1-based).

# Throws
- `BoundsError`: If `set_index` or `item_index` are out of valid bounds.
"""
function Base.setindex!(ss::SetSystem, value::Bool, set_index::Int, item_index::Int)
    if !(1 <= set_index <= ss.m_sets)
        error("Set index must be between 1 and $(ss.m_sets).")
    end
    if !(1 <= item_index <= ss.n_items)
        error("Item index must be between 1 and $(ss.n_items).")
    end
    
    # Use the item index as the row and set index as the column
    ss.sets[item_index, set_index] = value
end

function Base.getindex(ss::SetSystem, set_index::Int, item_index::Int)::Bool
    if !(1 <= set_index <= ss.m_sets)
        error("Set index must be between 1 and $(ss.m_sets).")
    end
    if !(1 <= item_index <= ss.n_items)
        error("Item index must be between 1 and $(ss.n_items).")
    end
    
    return  ss.sets[item_index, set_index]
end


"""
    greedy_set_cover(ss::SetSystem)

Implements the greedy set cover algorithm. It finds an append!roximate solution
to the set cover problem by iteratively selecting the set that covers the most
uncovered items.

# Arguments
- `ss::SetSystem`: The SetSystem instance for which to find the set cover.

# Returns
- `Vector{Int}`: A vector containing the indices of the selected sets
  that form the greedy set cover.
"""
function greedy_set_cover(ss::SetSystem)
    # The set of all items that need to be covered
    universe = Set(1:ss.n_items)
    
    # The set of items that have been covered so far
    covered_items = Set{Int}()
    
    # The list of sets selected to be in the cover
    selected_sets = Int[]
    
    
    # A boolean array to track which sets are still available
    sets_available = fill(true, ss.m_sets)

    # Continue until all items in the universe are covered
    while length(covered_items) < ss.n_items
        best_set_index = -1
        max_uncovered_count = 0
        
        # Iterate through all sets to find the one that covers the most new items
        for j in 1:ss.m_sets
            if sets_available[j]
                current_uncovered_count = 0
                # Iterate through all items
                for i in 1:ss.n_items
                    # If the item is in the current set and not yet covered
                    if ss.sets[i, j] && !(i ∈ covered_items)
                        current_uncovered_count += 1
                    end
                end
                
                if current_uncovered_count > max_uncovered_count
                    max_uncovered_count = current_uncovered_count
                    best_set_index = j
                end
            end
        end

        # If no new items can be covered, break the loop
        if max_uncovered_count == 0
            break
        end

        # Add the best set to the solution
        push!(selected_sets, best_set_index)

        # Update the covered items with the new items from the selected set
        for i in 1:ss.n_items
            if ss.sets[i, best_set_index]
                push!(covered_items, i)
            end
        end

        # Mark the selected set as unavailable
        sets_available[best_set_index] = false
    end
    
    return selected_sets
end


####################################################################################

function  compute_all_anchor_boxes( PS::Vector{Point2F}, ε )
    return  [sq_from_anchor( p, ε ) for p ∈ PS]
end

function  compute_set_system( PS::Vector{Point2F}, B::Vector{BBox2F} )
    ss = SetSystem( length( PS ), length( B ) )
    for i ∈ 1:length( PS )
        for  set ∈ 1:length( B )
            if  ( PS[ i ] ∈ B[ set ] )
                ss[ set, i ] = true
            end
        end
    end

    return  ss
end


function  greedy_cover( PS::Vector{Point2F}, ε )
    B = compute_all_anchor_boxes( PS, ε )
    ss = compute_set_system( PS, B )
    sets = greedy_set_cover( ss )

    return  [B[ s ]  for  s ∈ sets]
end


function  cleanup_cover( PS, cover )
    ss = compute_set_system( PS, cover )
    sets = greedy_set_cover( ss )

    return  [cover[ s ]  for  s ∈ sets]
end


####################################################################################

function  create_cover_3_approx( _PS::Vector{Point2F}, ε )
    PS = deepcopy( _PS )

    sort!( PS, by=p -> (p[1],-p[2]), rev=true )

    boxes = Vector{BBox2F}()
    alive = fill(true, length( PS ) )
    for  i ∈ 1:length( PS )
        if  !alive[ i ] 
            continue
        end
        p = PS[ i ]
        new_boxes = cover_for_point( p, ε )
        for  bb ∈ new_boxes
            f_add = false
            for j ∈ 1:length(PS)
                if  alive[ j ]  &&  ( PS[ j ] ∈ bb )
                    alive[ j ] = false
                    f_add = true
                end
            end
            if  ( f_add )
                push!( boxes, bb )
            end
        end
    end
    return  boxes
end


# Generate the upper grid 
function  lift_to_2d( P::Polygon1F )
    PS = Vector{Point2F}()
    n = length( P )
    for  i ∈ 1:n
        for  j ∈ (i+1):n
            push!( PS, Point2F( P[i][1], P[j][1] ) )
        end
    end
    return  PS
end
        


function  WSPD_cover( P::Polygon1F, ε )
    W = WSPD.init( P, ε )
    WSPD.expand( W )

    boxes = Vector{BBox2F}()
    count = 0
    for  pair ∈ W.finals
        L, R = WSPD.get_orig_ranges( W, pair )
        if  ( first( L ) > first( R ) )
            L,R = R,L
        end
        if  ! ( first( L ) <= last( L ) < first( R) <= last( R ) )
            println( L )
            println( R )
            
            @assert( first( L ) <= last( L ) < first( R) <= last( R ) )
        end
        γ = 0.0001
        p_l = Point2F( first( L )-γ, first( R )-γ )
        p_r = Point2F( last( L ) + γ, last( R ) + γ )
        bb = BBox2F_init( p_l, p_r )
        #bbox.expand( bb, 0.0001 )

        count += 1
        println( count, " : " )
        println( bb )
        push!( boxes, bb )
    end
    return  boxes
end


function     draw_cover( cr, PS, _boxes, radius = 0.15, f_expand = true )    
    # Expand each box in _boxes by 0.3 in all directions...
    B = Vector{BBox2F}()
    for  bb ∈ _boxes
        if  f_expand 
            if  cg.width( bb ) < 0.1  &&  cg.height( bb ) < 0.1
                push!( B, bb + 0.3 )
            else
                push!( B, bb + 0.1 )
            end
        else
            push!( B, bb )
        end
    end
    
    set_source_rgba( cr, 0.5, 1.0, 0.5, 0.5 )
    fill_bboxes_cycle_colors( cr, B, α=0.3 )
    
    set_source_rgba( cr, 0.0, 0.5, 0.0, 1.0 )
    draw_bboxes( cr, B ) 

    set_source_rgb( cr, get_color_rgb( 2 )... )
    set_line_width(cr, 1.0)
    draw_points( cr, PS, radius )
end


function     draw_union( cr, PS, _boxes, radius = 0.15, f_expand = true )    
    B = sort( _boxes, by=b -> energy( b ) )
    
    set_source_rgba( cr, 0.5, 1.0, 0.5, 1.0 )
    fill_bboxes_cycle_colors( cr, B; α=1.0, f_draw_frames = true )
    
    set_source_rgb( cr, get_color_rgb( 2 )... )
    set_line_width(cr, 1.0)
    draw_points( cr, PS, radius )
end


function compute_opt_cover_ip( PS, ε )
    B = compute_all_anchor_boxes( PS, ε )
    ss = compute_set_system( PS, B )

    items = length( PS )
    sets = length( B )

    model = Model(Gurobi.Optimizer)
    #model = Model(HiGHS.Optimizer)
    #model = Model(HiGHS.Optimizer)

    # Verify set system
    for e in 1:items
        for s in 1:sets
            if  ss[ s, e ]
                @assert( PS[ e ] ∈ B[ s ] )
            else
                @assert( ! (PS[ e ] ∈ B[ s ]  ) )
            end
        end
    end
    #exit( -1 )
    
    # The decision variables are binary (0 or 1), representing whether each subset
    # is chosen for the cover.
    @variable(model, x[1:sets], Bin)

    # The objective function is to minimize the total number of subsets chosen.
    @objective(model, Min, sum(x[s] for s in 1:sets ) )

    # The constraints ensure that every element in the universe is covered
    # by at least one selected subset.
    for  e ∈ 1:items
        @constraint(model, sum(x[s] for s in 1:sets if ss[ s, e ] ) >= 1
        )
    end
    #println( model )
    
    # Solve the optimization problem.
    optimize!(model)

    # Check the status of the solution.
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found.")

        # Get the value of the objective function (minimum number of sets).
        println("Minimum number of sets required: ", round(Int, objective_value(model)))
        
        # Get the values of the decision variables to see which sets were chosen.
        optimal_sets = [s for s in 1:sets if value(x[s]) > 0.5]
        println("Sets in the optimal cover: ", optimal_sets)
    else
        println("No optimal solution found. Termination status: ", termination_status(model))
    end
    return  [B[s] for s ∈ optimal_sets ]
end


function  split_by_dim( b::BBox{D,T}, dim::Int, val::T ) where{D,T}
    if  ( val < d_min( b, dim ) )  ||  ( val > d_max( b, dim ) )
        return  [b]
    end

    maxi = deepcopy( b.maxi )

    out = Vector{BBox{D,T}}()

    maxi[ dim ] = val
    push!( out, BBox{D,T}( Point{D,T}( b.mini... ), Point{D,T}( maxi... ) ) )
    
    mini = deepcopy( b.mini )
    mini[ dim ] = val
    push!( out, BBox{D,T}( Point{D,T}( mini... ), Point{D,T}( b.maxi... ) ) )

    return  out
end


#=
function  split_vertically( rect::BBox2F, x )
    x_range = get_dim( rect, 1 )
    y_range = get_dim( rect, 2 )

    b_1 = Bbox2F( Point2F( x_range[1], y_range[1] ), Point2F( x, y_range[2] ) ) 
    b_2 = Bbox2F( Point2F( x, y_range[1] ), Point2F( x_range[2], y_range[2] ) ) 

    return  b_1, b_2
end

#function  base.(×)( )
function  split_horizontally( rect::BBox2F, y )
    x_range = get_dim( rect, 1 )
    y_range = get_dim( rect, 2 )

    b_1 = Bbox2F( Point2F( x_range[1], y_range[1] ), Point2F( x_range[2], y ) ) 
    b_2 = Bbox2F( Point2F( x, y_range[1] ), Point2F( x_range[2], y_range[2] ) ) 

    return  b_1, b_2
end
=#
function  split_box( curr::BBox2F, p::Point2F, orect::BBox2F )
    sp = split_by_dim( curr, 1, p[1 ] )
    B = mapreduce( b -> split_by_dim( b, 2, p[2] ), vcat, sp )

    #=for  b ∈ B
        println( "b:", b,  "  ", b ⊆ orect )
    end=#
    u = filter( b -> ! ( b ⊆ orect ), B )

    #=println( length( u ) )
    println( "Stop for now" )
    exit( -1 )=#
    return  u
end


function  extract_ranges( V::Vector{Point{D,T}}, s, dim ) where {D,T} 
    # x axes needed splitting
    C = unique!(sort( [v[ dim ] for v ∈ V] ) )
    X = filter( x -> ( d_min( s, dim ) < x < d_max(s,dim ) ), C ) 
    
    pushfirst!(X, d_min( s, dim ) )
    push!(X, d_max( s, dim ) )
    unique!( sort!( X ) )

    out = Vector{Tuple{T,T}}()
    for  i ∈ 1:(length(X) - 1)
        push!( out, (X[i],X[i+1]) )
    end
    return  out
end


function  is_covered( bb::BBox2F, S::Vector{BBox2F} )
    for  s ∈ S
        if  bb ⊆ s
            return  true
        end
    end
    return  false
end

function  print_rects( O )
    println( "\n/-------------------------" );
    for r ∈ O
        println( r )
    end
    println( "\\-------------------------\n" );
end



"""

splinter the square s into rectangles, by the vertices of the squares
of S that lies inside s.

"""
function  splinter( s::BBox2F, S::Vector{BBox2F} )
    V = Vector{Point2F}() # All vertices of S inside s
    for  curr ∈ S
        if  is_intersect( curr, s )
            Base.append!( V, vertices( curr ) )
        end
    end
    if length( V) == 0
        return  [s]
    end

    X =  extract_ranges( V, s, 1 )
    Y =  extract_ranges( V, s, 2 )

    O = Vector{BBox2F}()
    for  x_r ∈ X
        sout = Vector{BBox2F}()
        for  y_r ∈ Y
            bb = x_r * y_r
            if  is_covered( bb, S )  ||  isempty( bb )
                continue
            end
            if  !isempty( sout )
                f_com, _ = compatible( last( sout ), bb )
                if  f_com
                    bbo = pop!( sout )
                    bb = merge( bb, bbo )
                end
            end
            push!( sout, bb )
        end
        if  (! isempty( sout ))
            Base.append!( O, sout )
        end
    end

    if  length( O ) <= 1
        return  O
    end
    
    ### We need to merge rectangles that are compatible in the y axis...
    sort!( O, by=r -> ( d_min( r, 2 ), d_max( r, 2 ) ) )
    f_found = false
    for  i ∈ 1:length( O ) - 1 
        for  j ∈ i+1:length( O ) 
            f_comp, _ = compatible( O[i], O[ j ] )
            if  f_comp
                f_found = true
            end
        end
    end

    
    OO = Vector{BBox2F}()
    f_f_found = false;
    f_comp::Bool = false
    old = O[1]
    for  i ∈ 2:length( O )
        f_comp, _ = compatible( old, O[i] )
        if  !f_comp
            push!(OO, old )
            old = O[ i ]
        else
            old = merge( old, O[i ] )
            f_f_found = true
        end
    end
    push!(OO, old )
    @assert( f_found == f_f_found )

    for  b ∈ OO
        @assert( b ⊆ s )
    end
    
    return  OO
end


"""
compute_disjoint_union

Decompose a set of shadow squares into disjoint rectangles. It is
*important* that these are shadow squares, as the logic works only for
them.
"""
function   compute_disjoint_union( shadow_squares )
    B = sort( shadow_squares, by=b -> -energy( b ) )
    rects_out = Vector{BBox2F}()
    handled = Vector{BBox2F}()

    count = 0;
    for  b ∈ B
        count += 1
        print( count, "     \r" )
        flush( stdout )
        sp = splinter( b, handled )
        Base.append!( rects_out, sp )
        push!( handled, b )
    end

    return  rects_out
end


function  random_anchored_squares( n, N, ε )
    B = Vector{BBox2F}()
    AN = Vector{Point2F}()
    for i ∈ 1:N
        anchor = Point2F( rand() * n,  rand() * n )
        if  anchor[ 2 ] < anchor[1]
            anchor = Point2F( anchor[2], anchor[1] )
        end
        push!( AN, anchor )
        push!( B, sq_from_anchor( anchor, ε ) )
    end
    return  B, AN
end

function one_dim_comp_solution( ε = 1.0, n = 80 )
    P = Polygon1F()
    for i ∈ 1:n
        push!( P, Point1F( i ) )
    end
    
    # Generate the upper grid 
    PS = lift_to_2d( P )

    cover_wspd = WSPD_cover( P, ε )
    cover_3_aprx = create_cover_3_approx( PS, ε )
    cover_greedy = greedy_cover( PS, ε )

    cover_3_aprx_c = cleanup_cover( PS, cover_3_aprx )

    println( "Running IP solver..." )
    cover_ip = compute_opt_cover_ip( PS, ε ) 

    N = length( PS )

    println( "n = ", n, " N = ", N )
    println( "WSPD        : ", length( cover_wspd ) )
    println( "3-approx    : ", length( cover_3_aprx ) )
    println( "3-approx(c) : ", length( cover_3_aprx_c ) )
    println( "Greedy      : ", length( cover_greedy ) )
    println( "IP          : ", length( cover_ip ) )

    #####################################################################
    # Drawing stuff... 
    c,cr,bb_draw = cairo_setup( "out/1_dim.pdf", Point2F( n+1, n+1 ) )
    flip_y_axis( cr, n+1 )

    str = @sprintf( "WSPD cover\nn = %d\nε = %g\nN = %d  (# of grid points)\n|WSPD| = %d\n",
                      n, ε, N, length( cover_wspd ) )
    add_to_pdf_text( cr, n, str )
    Cairo.show_page( cr )
    
    draw_cover( cr, PS, cover_wspd )
    Cairo.show_page( cr )

    ##########################################################
    str = @sprintf( "3 approx cover\nn = %d\n|3-approx| = %d\n",
                      n,  length( cover_3_aprx ) )
    add_to_pdf_text( cr, n, str )

    draw_cover( cr, PS, cover_3_aprx )
    Cairo.show_page( cr )

    ##########################################################
    str = @sprintf( "3 approx cover (clean)\nn = %d\n|3-approx-clean| = %d\n%s",
        n,  length( cover_3_aprx_c ),
        "Cleaning is done using greedy cover on the 3-approx\n"*
        "solution computed." )
    add_to_pdf_text( cr, n, str )
    Cairo.show_page( cr )

    draw_cover( cr, PS, cover_3_aprx_c )
    Cairo.show_page( cr )

    ##########################################################
    str = @sprintf( "Greedy cover\nn = %d\nε = %g\n|Greedy| = %d\n%s",
        n, ε, length( cover_greedy ),
        "Solution computed by setting up the set system,\n"*
        "and using greedy set-cover to solve it." )
    add_to_pdf_text( cr, n, str )
    Cairo.show_page( cr )

    draw_cover( cr, PS, cover_greedy )
    Cairo.show_page( cr )

    ##########################################################
    str = @sprintf( "IP cover\nn = %d\nε = %g\n|IP| = %d\n%s",
        n, ε, length( cover_ip ),
        "Set up the natural IP for set-cover, and solve it exactly\n"
        *"using Gurobi" )
    add_to_pdf_text( cr, n, str )
    Cairo.show_page( cr )

    draw_cover( cr, PS, cover_ip )
    Cairo.show_page( cr )

    ##########################################################
    str = @sprintf( "Summary\nn = %d    (N=%d)\nε = %g\n",
        n, N, ε )
    str_ip       = @sprintf( "|IP|       = %d\n", length( cover_ip ) )
    str_greedy   = @sprintf( "|Greedy|   = %d\n", length( cover_greedy ) )
    str_wspd     = @sprintf( "|WSPD|     = %d\n", length( cover_wspd ) )
    str_r_aprx   = @sprintf( "|3-aprx|   = %d\n", length( cover_3_aprx ) )
    str_r_aprx_c = @sprintf( "|3-aprx-c| = %d\n", length( cover_3_aprx_c ) )

    sout = str * str_ip * str_greedy * str_wspd * str_r_aprx * str_r_aprx_c
    add_to_pdf_text( cr, n, sout )

    ##########################################################
    Cairo.finish( c )
end

function  draw_diagonal( cr, n )
    set_line_width(cr, 1.0)
    Cairo.move_to(cr, 0, 0 )
    Cairo.line_to(cr, n+1, n+1 )
    Cairo.stroke( cr )
end

function one_dim_union(; ε = 1.0, n = 80 )
    P = Polygon1F()
    for i ∈ 1:n
        push!( P, Point1F( i ) )
    end
    
    # Generate the upper grid 
    PS = lift_to_2d( P )
    
    cover_greedy = greedy_cover( PS, ε )
    cover_wspd = WSPD_cover( P, ε )

    println( "COVERS computed" )
    N = length( PS )

    println( "n = ", n, " N = ", N )
    println( "Greedy      : ", length( cover_greedy ) )

    #####################################################################
    # Drawing stuff... 
    fln = "out/union.pdf"
    c,cr,_ = cairo_setup( fln, Point2F( n+1, n+1 ) )
    flip_y_axis( cr, n+1 )

    # Draw the greedy solution
    str = @sprintf( "%sn = %d\nε = %g\nN = %d  (# of grid points)\nsize = %d\n",
                    "Greedy cover\n",
                      n, ε, N, length( cover_greedy ) )
    add_to_pdf_text( cr, n, str )

    draw_cover( cr, PS, cover_greedy, 0.1, false )    

    Cairo.show_page( cr )

    println( "A" );
    ####################################################################
    str = @sprintf( "%sn = %d\nε = %g\nN = %d  (# of grid points)\nsize = %d\n",
                    "Greedy cover drawn from smallest to largest square\n",
                      n, ε, N, length( cover_greedy ) )
    add_to_pdf_text( cr, n, str )

    draw_union( cr, Vector{Point2F}(), cover_greedy )
    
    Cairo.show_page( cr )
    println( "B" );

    # Now compute the union
    rects_greedy = compute_disjoint_union( cover_greedy )
    str = @sprintf( "%s%d ==> %d   (WSPD solution size: %d\n",
                    "Greedy solution broken into disjoint rectangles\n",
        length( cover_greedy),  length( rects_greedy  ),
        length( cover_wspd )
    )
    add_to_pdf_text( cr, n, str )

    println( "C" );
    #draw_cover( cr, Vector{Point2F}(), rects_greedy, 0.0, false )
    draw_cover( cr, PS, rects_greedy, 0.1, false )
    Cairo.show_page( cr )
    

    
    ###############################################################
    # Now for a more interesting set of squares...
    sz = div( n * n, 4 )
    B, AN = random_anchored_squares( n, sz , ε )
    @assert( length( B ) ==  sz )
    str = @sprintf( "%s\n|size| = %d\n",
                    "Random collection of shadowed squares\n",
                      length( B ) )
    add_to_pdf_text( cr, n, str )

    println( "D" );
    draw_union( cr, Vector{Point2F}(), B )
    draw_diagonal( cr, n )

    set_source_rgb( cr, get_color_rgb( 2 )... )
    set_line_width(cr, 1.0)
    draw_points( cr, AN, 0.07 )
    
    println( "E" );
    Cairo.show_page( cr )
    ###############################################################

    println( "G:", length( B ) );
    rects = compute_disjoint_union( B )
    str = @sprintf( "%s\n %d ==>  = %d\n",
                    "Random collection of shadowed squares as disjoint union\n",
                     length( B), length( rects ) )
    add_to_pdf_text( cr, n, str )
    println( "F" );

    draw_cover( cr, Vector{Point2F}(), rects, 0.0, false )
    
    Cairo.finish( c )

    println( "Created : ", fln )
end

#########################################################################



    
function  draw_vicinity_point( cr, p, ε )
    PS = Vector{Point2F}()
    push!( PS, p )
    
    sq_p = sq_from_anchor( p, ε )
    sq_p_tl = sq_from_top_left( p, ε )
    sq_p_tr = sq_from_top_right( p, ε )
    sq_p_bl = sq_from_bottom_left( p, ε )
    
    B = [sq_p sq_p_tl sq_p_tr sq_p_bl]
    
    # A bit of Julia cleverness
    V = vcat( ( vertices(b) for b ∈ B )... )
    
    P = Polygon2F( convex_hull_2d( V ) )

    #Cairo.set_source_rgb(cr, 1.0, 0.0, 0.0) # RGB for red
    set_source_rgba( cr, 1.0, 0.0, 0.0, 0.1 )
    fill_polygon( cr, P )
    
    set_source_rgb( cr, get_color_rgb( 1 )... )
    draw_polygon( cr, P, true )
        
    draw_cover( cr, PS, B, radius = 0.05, f_expand = false ) #_boxes, f_expand = true )    
end


function   draw_vicinity( filename = "out/vicinity.pdf", ε = 1.0 )
    n = 10    
    fln = "out/vicinity.pdf"
    c,cr,bb_draw = cairo_setup( fln, Point2F( n+1, n+1 ) )

    flip_y_axis( cr, n+1 )

    set_source_rgb( cr, get_color_rgb( 2 )... )
    set_line_width(cr, 1.0)
    Cairo.move_to(cr, 0, 0 )
    Cairo.line_to(cr, n+1, n+1 )
    Cairo.stroke( cr )

    draw_vicinity_point( cr, Point2F( 3.0, 5.0 ), ε )
    draw_vicinity_point( cr, Point2F( 3.0, 8.5 ), ε )
    draw_vicinity_point( cr, Point2F( 5.0, 6.0 ), ε )
    draw_vicinity_point( cr, Point2F( 8.0, 9.0 ), ε )
    draw_vicinity_point( cr, Point2F( 9.0, 10.0 ), ε )
    Cairo.show_page( cr )

    p = Point2F(4.0, 8.5 )
    draw_vicinity_point( cr, p, ε )
    sq_p_tr = sq_from_top_right( p, ε )
    q = bottom_left( sq_p_tr )
    sq_q = sq_from_anchor( q, ε )
    
    set_source_rgba( cr, 0.5, 0.5, 0.0, 0.3 )
    fill_bboxes( cr, [sq_q] ) 
    set_source_rgba( cr, 0.0, 0.0, 0.0, 1.0 )
    draw_bboxes( cr, [sq_q] ) 
    set_source_rgba( cr, 0.0, 0.0, 1.0, 1.0 )
    draw_points( cr, [q], 0.05 )
    Cairo.show_page( cr )
    
    Cairo.finish( c )

    println( "Created : ", fln ) 
end

function (@main)(ARGS)    
    if length( ARGS ) > 0  &&  ARGS[1] == "1dim"
        one_dim_comp_solution( ε=0.9, n = 20 )
        return
    end

    if  ( length( ARGS ) > 0  &&  ARGS[1] == "union" )
        one_dim_union( ε=0.7123, n = 80 )
        return
    end

    if length( ARGS ) > 0  &&  ARGS[1] == "vicinity"
        draw_vicinity( "vicinity.pdf", ε = 0.4 )
    end
    return  0
end


####################################################################


#! /bin/env  julial

# shortcut.jl:
#
# Approximate the optimal shortcut for a give data file.
#

push!(LOAD_PATH, pwd()*"/src/")
#push!(LOAD_PATH, pwd()*"/src/cg/")

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.bbox;
using FrechetDist.cg.point;
using VirtArray;
using Printf;
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

Point1F = Point{1,Float64};
Polygon1F = Polygon{1,Float64};

########################################################################################
########################################################################################

# The anchor stuff

function  energy( p::Point2F )::Float64
    @assert( p[2] >= p[ 1 ] )
    return  p[ 2 ] - p[ 1 ]
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
    ℓ = ( energy( p ) - new_e ) / 2.0;

    new_anc = Point2F( x + ℓ, y - ℓ )
    return  sq_from_anchor( new_anc, ε )
end

function  sq_from_bottom_left( p::Point2F, ε )
    x, y = p[ 1 ], p[2] 
    new_e = energy( p )  / ( 1.0 + ε ) 
    ℓ = energy( p ) - new_e;

    new_anc = Point2F( x + ℓ, y )
    return  sq_from_anchor( new_anc, ε )
end


function  sq_from_top_right( p::Point2F, ε )
    x, y = p[ 1 ], p[2] 
    new_e = energy( p )  / ( 1.0 + ε ) 
    ℓ = energy( p ) - new_e;

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

Implements the greedy set cover algorithm. It finds an approximate solution
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

    return  ss;
end


function  greedy_cover( PS::Vector{Point2F}, ε )
    B = compute_all_anchor_boxes( PS, ε )
    ss = compute_set_system( PS, B )
    sets = greedy_set_cover( ss )

    return  [B[ s ]  for  s ∈ sets];
end


function  cleanup_cover( PS, cover )
    ss = compute_set_system( PS, cover )
    sets = greedy_set_cover( ss )

    return  [cover[ s ]  for  s ∈ sets];
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
        p = PS[ i ];
        new_boxes = cover_for_point( p, ε )
        for  bb ∈ new_boxes
            f_add = false;
            for j ∈ 1:length(PS)
                if  alive[ j ]  &&  ( PS[ j ] ∈ bb )
                    alive[ j ] = false;
                    f_add = true;
                end
            end
            if  ( f_add )
                push!( boxes, bb )
            end
        end
    end
    return  boxes;
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
    count = 0;
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
    return  boxes;
end


function     draw_cover( cr, PS, _boxes; radius = 0.15, f_expand = true )    
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
    fill_bboxes_cycle_colors( cr, B, 0.3 )
    
    set_source_rgba( cr, 0.0, 0.5, 0.0, 1.0 )
    draw_bboxes( cr, B ) 

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
    return  [B[s] for s ∈ optimal_sets ];
end


function add_to_pdf_text( cr, n, text_to_write::String )
    Cairo.save( cr )
    flip_y_axis( cr, n+1 )

    scale = n * 0.02 / 8.0;
    Cairo.scale(cr, scale, scale )
    set_source_rgba( cr, 0.0, 0.0, 0.0, 1.0 )

    
    # Set the font properties
    font_size = 14.0
    Cairo.set_font_size(cr, font_size)

    # Set the font face. You must have this font installed on your system.
    # We will use "Serif" as a common fallback.
    Cairo.select_font_face(cr, "Serif", Cairo.FONT_SLANT_NORMAL,
        Cairo.FONT_WEIGHT_NORMAL)

    # Split the input string by newline characters
    lines = split(text_to_write, '\n')

    # Define the starting position and line spacing
    x_pos, y_pos = 10.0, 10.0
    line_spacing = font_size * 1.5

    # Iterate through each line and draw it
    for line in lines
        # Move to the position for the current line
        Cairo.move_to(cr, x_pos, y_pos)
        # Show the text
        Cairo.show_text(cr, line)
        # Increment the vertical position for the next line
        y_pos += line_spacing
    end
    Cairo.show_page( cr )
    Cairo.restore( cr )
end


function one_dim_comp_solution(; ε = 1.0, n = 80 )
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

    draw_cover( cr, PS, cover_3_aprx_c )
    Cairo.show_page( cr )

    ##########################################################
    str = @sprintf( "Greedy cover\nn = %d\nε = %g\n|Greedy| = %d\n%s",
        n, ε, length( cover_greedy ),
        "Solution computed by setting up the set system,\n"*
        "and using greedy set-cover to solve it." )
    add_to_pdf_text( cr, n, str )

    draw_cover( cr, PS, cover_greedy )
    Cairo.show_page( cr )

    ##########################################################
    str = @sprintf( "IP cover\nn = %d\nε = %g\n|IP| = %d\n%s",
        n, ε, length( cover_ip ),
        "Set up the natural IP for set-cover, and solve it exactly\n"
        *"using Gurobi" )
    add_to_pdf_text( cr, n, str )

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

    sout = str * str_ip * str_greedy * str_wspd * str_r_aprx * str_r_aprx_c;
    add_to_pdf_text( cr, n, sout )

    ##########################################################
    Cairo.finish( c )


    return  0;
end

function  vertices( bb::BBox2F )
    return  [ bottom_left( bb ), bottom_right( bb ),
              top_left( bb ), top_right( bb ) ]
end


    
function  draw_vicinity_point( cr, p, ε )
    PS = Vector{Point2F}()
    push!( PS, p )
    
    sq_p = sq_from_anchor( p, ε )
    sq_p_tl = sq_from_top_left( p, ε )
    sq_p_tr = sq_from_top_right( p, ε )
    sq_p_bl = sq_from_bottom_left( p, ε )
    
    B = [sq_p sq_p_tl sq_p_tr sq_p_bl];
    
    # A bit of Julia cleverness
    V = vcat( ( vertices(b) for b ∈ B )... )
    
    P = Polygon2F( convex_hull_2d( V ) )

    #Cairo.set_source_rgb(cr, 1.0, 0.0, 0.0) # RGB for red
    set_source_rgba( cr, 1.0, 0.0, 0.0, 0.1 )
    fill_polygon( cr, P )
    
    set_source_rgb( cr, get_color_rgb( 1 )... )
    draw_polygon( cr, P, true );
        
    draw_cover( cr, PS, B, radius = 0.05, f_expand = false ) #_boxes, f_expand = true )    
end


function   draw_vicinity( filename = "vicinity.pdf"; ε = 1.0 )
    n = 10    
    c,cr,bb_draw = cairo_setup( "out/vicinity.pdf", Point2F( n+1, n+1 ) )

    flip_y_axis( cr, n+1 )

    set_source_rgb( cr, get_color_rgb( 2 )... )
    set_line_width(cr, 1.0)
    Cairo.move_to(cr, 0, 0 )
    Cairo.line_to(cr, n+1, n+1 )
    Cairo.stroke( cr );

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
end


function (@main)(ARGS)
    if length( ARGS ) > 0  &&  ARGS[1] == "1dim"
        one_dim_comp_solution( ε=0.9, n = 20 )
        return
    end

    draw_vicinity( "vicinity.pdf", ε = 0.4 )
    return  0;
end

####################################################################3

#export  BBTree_build;

#end

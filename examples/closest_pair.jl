#! /bin/env julia

########################################################################
# Closest pair of n points in constant dim in linear time using
# hashing, based on Rabin's algorithm (simplification of this
# algorithms to be precise).
#
# Implemented by Sariel Har-Peled
########################################################################
# I invested some energy in making the program being more
# efficient. One heuristic that the grid is recomputed only if the
# min-distance changed, and the load of some cell exceeds 4
#
# A minor technical issue is that Julia does not support "open"
# hashing directly. To overcome this, every grid cell has a vector of
# all the ids of the points stored in this cell. This seems wasteful,
# but the performance is decent enough. I am sure a speedup of 30%-50%
# might be possible by doing this in a more clever way, but the
# solutions I tried (and deleted he code from here), did no in fact
# yield significant speedup considering their complexity.
#
# things in Julia are a bit weird. The exact way you fill the array of
# points by random points can have a speed difference of up to
# 10x. The current implementation (essentially using the library
# function correct, seems to be fastest.
#
# 2025-June-7
########################################################################


push!(LOAD_PATH, pwd()*"/src/")
push!(LOAD_PATH, pwd()*"/src/cg/")

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using Printf;
using StaticArrays;

######################################################################
# Misc
######################################################################

@inline function  gid(P::Point{D,T}, r::T)::Point{D,Int} where{D,T}
    x = MVector{D,Int}(undef);

    for i ∈ 1:D
        x[i] = floor( Int, P[ i ] / r );
    end

    return  Point{D,Int}( x );
end


######################################################################
# Linear time algorithm using hashing - using vectors to implement
# open hashing - probably a terrible idea.
######################################################################


mutable struct GridType{D,T}
    P::Vector{Point{D,T}}
    ℓ::Float64  # Side length
    cells::Dict{Point{D,Int},Vector{Int}}

    cp::Tuple{Int, Int};
    cp_dist::T;
    f_regrid::Bool

    one::Point{D,Int}  # A vector filled with ones, oh yeh!
end


"""
    add_value!

    An insertion function for dictionary where a value is an array of values.
"""
@inline function add_value!(dict::Dict{Point{D,Int},Vector{Int}},
                            key::SVector{D,Int}, value::Int) where {D}
    if haskey(dict, key)
        push!(dict[key], value)
    else
        dict[key] = [value]
    end
end

function   grid_init( P::Vector{Point{D,T}}, ℓ::T,
                      r::UnitRange{Int} = eachindex( arr ), cp = (-1,-1) ) where{D,T}

    dict = Dict{Point{D,Int},Vector{Int}}();
    sizehint!( dict, min( length( r ), length( P ) ) )
    G = GridType( P, ℓ, dict, cp, ℓ, false,
                  Point{D,Int}( fill( 1, D ) ), MVector{D,Int}(undef ) );

    for  i ∈ r
        trg = gid(P[i], ℓ )
        add_value!( G.cells, trg, i );
    end

    return  G;
end

function  closest_pair_add_point( G::GridType{D,T}, loc::Int ) where{D,T}
    p =  G.P[ loc ];
    trg = gid( p, G.ℓ )
    id_min = trg - G.one;
    id_max = trg + G.one;

    P = G.P;
    min_ind = -1;
    min_dist = G.cp_dist;
    for  cell ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        i = Point{D,Int}( Tuple( cell )... );
        #println( typeof( i ) )
        if  ! haskey( G.cells, i )  continue  end
        list = G.cells[ i ];
        if  length( list ) > 4
            G.f_regrid = true;
        end
        for  p_ind ∈ list
            new_dist = Dist( P[ p_ind ], p )
            if   new_dist < min_dist
                min_dist = new_dist;
                min_ind = p_ind;
            end
        end
    end

    # Minimum distance had not changed. Nothing much to do...
    add_value!( G.cells, gid( P[ loc ], G.ℓ ), loc );
    if   min_ind > 0
        G.cp = (min_ind, loc );
        G.cp_dist = min_dist;
    end

    if  G.f_regrid  
        #println( "REGRID!" );
        G = grid_init( G.P, G.cp_dist, 1:loc, G.cp );
    end

    return  G;
end


function  closest_pair( P::Vector{Point{D,T}} ) where{D,T}
    d = Dist( P[ 1 ], P[ 2 ] );
    @assert( d > 0.0 );

    G = grid_init( P, d, 1:2, (1,2) );
    for  i ∈ 3:length( P )
        G = closest_pair_add_point( G, i );
    end

    return  G;
end


############################################################################
# Closest pair in the plane using divide-and-conquer in O( n log n ) time.
############################################################################

mutable struct  Solution
    d::Float64
    i::Int64
    j::Int64
end
function  Sol_set( sol::Solution, _d, _i, _j )
    sol.d, sol.i, sol.j = _d, _i, _j;
end

Points2F = Vector{Point2F};

struct  PointSetXY
    P::Points2F;
    ord_x::Vector{Int64};
    ord_y::Vector{Int64};
end

len( PS::PointSetXY ) = length( PS.ord_x );    # number of points in PS
#d( p, q ) =  norm( ( p[1] - q[1], p[2] - q[2] ) ); # Distance between points
PS_ord_y( Q, l ) = (Q.P[ Q.ord_y[ l ] ])[ 2 ];

struct SortByCoord <: Base.Order.Ordering
    P::Vector{Point2F};
    coord::Int
end

import Base.Order.lt
function  lt(o::SortByCoord, a, b)
    p,q = o.P[ a ], o.P[ b ];
    return  (p[o.coord] < q[o.coord] );
end

function  PointSetXY( _P::Points2F )
    #PS = PointSetXY( _P, collect( 1:length( _P ) ), collect( 1:length( _P ) ) );
    len::Int = length( _P );
    v1n=[i for i ∈ 1:len ]
    PS = PointSetXY( _P, v1n, deepcopy( v1n ) );

    # Sort the points by x-axis, and store the ordering in PS.ord_x
    sort!( PS.ord_x, order=SortByCoord( PS.P, 1 ) );

    # Sort the points by y-axis, and store the ordering in PS.ord_y
    sort!( PS.ord_y, order=SortByCoord( PS.P, 2 ) );

    return  PS;
end

# Computing subset of points in _PS that fulfill the condition f.
function  PS_filter( P::PointSetXY, f )
    return    PointSetXY( P.P, filter( i -> f( P.P[ i ] ), P.ord_x ),
                             filter( i -> f( P.P[ i ] ), P.ord_y ) );
end

# Use the elevator algorithm to discover any closest pair that is
# better than the current solution.
function  nn_middle( VPS::PointSetXY, mid, sol::Solution )
    low, hi = mid - sol.d, mid + sol.d
    P = PS_filter( VPS, p -> ( low <= p[1] <= hi ) );
    ( len( P ) < 2 )  &&  return  sol;

    b = t = 1;  # Bottom/top of elevator
    while  ( t <= len( P ) )
        # Delete bottom point ∈ elevator if irrelevant
        while  ( b < t ) && ( PS_ord_y( P, b ) < (PS_ord_y( P, t ) - sol.d) )
            b = b + 1;
        end

        # Find nn in points in the elevator to the top point in it.
        # The elevator points are all points with indices in
        # PS.ord_y[ b:t ]
        i_p = P.ord_y[ t ];
        p = P.P[ i_p ];

        for  j ∈ b:(t-1)
            i_q = P.ord_y[ j ]
            q = P.P[ i_q ]
            ℓ = Dist( p, q )
            if  ( ℓ < sol.d )
#                println( "XXXYY" );
                Sol_set( sol, ℓ, i_p, i_q );
            end
        end
        t = t + 1;
    end

    return  sol;
end

function  closest_pair_dc_inner( PS::PointSetXY, sol )
    ( len( PS ) < 2)  &&  return  sol;

    mid = PS.ord_x[ len( PS )  ÷  2 ];
    x_mid = PS.P[ mid ][ 1 ];

    # A subtlety: All the points with x-coordinat x_mid, would be
    # checked by the call nn_middle, so no need to include them in the
    # two recursive calls.

    # Dividing...
    PS_L = PS_filter( PS, p -> ( p[1] < x_mid ) );  # Left points
    PS_R = PS_filter( PS, p -> ( p[1] > x_mid ) );  # Right points

    # Recursing
    closest_pair_dc_inner( PS_L, sol ); # recurse left
    closest_pair_dc_inner( PS_R, sol ); # recurse right

    # Conquering
    nn_middle( PS, x_mid, sol );

    # declaring victory, and withdrawing
    return  sol;
end

function  closest_pair_dc( _PS::Points2F )
    PS = PointSetXY( _PS )

    sol = Solution( Dist(PS.P[1], PS.P[2]), 1, 2 );
    return  closest_pair_dc_inner( PS, sol );
end


############################################################################
# Closest pair  brute force in O(n^2) time.
############################################################################

function  closest_pair_brute_force( PS::Vector{Point{D,T}} ) where{D,T}
    sol = Solution( Dist(PS[1], PS[2]), 1, 2 );

    for  i in 1:length( PS ) -1
        p = PS[ i ];
        for  j in i+1:length( PS )
            q = PS[ j ];
            ℓ =  Dist( p, q )
            if  ( ℓ < sol.d )
                Sol_set( sol,  ℓ, i, j );
            end
        end
    end

    return  sol;
end


function  append_to_file( file_path, str )
    try
        open(file_path, "a") do file
            write(file,  str)
        end
    catch e
        println("An error occurred while appending to the file: ", e)
    end
end

############################################################################
############################################################################


function  force_compile()
    D=2
    n = 1000
    P = Points( Polygon_random( D,Float64, n ) );
    t_dc      = @timed sol = closest_pair_dc( P );
    t_bf      = @timed sol_bf = closest_pair_brute_force( P );
    t_rand    = @timed G = closest_pair( P );
end


function (@main)(ARGS)
    if  length( ARGS ) != 1
        println( "Usage:\n\t"
                 * "closest_pair.jl [n]\n\n" );
        return  -1;
    end

    D = 2;
    n =  parse(Int64, replace( ARGS[ 1  ], "," => "" ) );

    println( "n: ", n );
    force_compile();

    println( "\n\n--------------------\n\n" );
    println( "Generating input & Copying..." );
    P = Polygon_random( D,Float64, n );

    PB = deepcopy( Points( P ) );
    PC = deepcopy( Points( P ) );
    PD = deepcopy( Points( P ) );

    println( "Done." );

    f_brute_force::Bool = false;
    if  ( n < 10_000 )
        f_brute_force = true;
        t_brute = @timed sol_bf = closest_pair_brute_force( PB );
        @printf( "Time brute         : %10.5f\n", t_brute.time );
    end
    
    t_rand    = @timed G      = closest_pair( PC );
    @printf( "Time Rand          : %10.5f\n", t_rand.time );

    t_dc      = @timed sol    = closest_pair_dc( PD );
    @printf( "Time D&C           : %10.5f\n", t_dc.time );

    println( "Distance                            : ", sol.d );

    str = @sprintf( "%12d,  %11.6f, %11.6f\n",
        length( P ),
        t_dc.time, t_rand.time );

    append_to_file( "results/closest_pair.txt", str );

    if  ( sol.d != G.cp_dist )
        println( "\n\n\nBUG!!!!\n\n\n\n" );a
    end

    return  0;
end

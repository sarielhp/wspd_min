#! /bin/env julia


########################################################################
# Computes a packing (i.e., a net) of n points in constant dim in
# linear time using hashing, based on Har-Peled and Raichel Net &
# Prune paper.
#
# Implemented by Sariel Har-Peled
########################################################################
#
# 2025-June-7
########################################################################


myBase = "/home/sariel/prog/25/wspd"

push!(LOAD_PATH, myBase * "/src/")
push!(LOAD_PATH, myBase * "/src/cg/")

### To make emacs lsp mode works correctly...
macro ignore(args...) end

@ignore    include("../src/cg/FrechetDist.jl")
@ignore    include("../src/cg/polygon.jl")
@ignore    include("../src/cg/point.jl")
@ignore    include("../src/cg/cg.jl")

using Distributions
using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using Printf;
using StaticArrays;
using Cairo

include( "graphics.jl" )

@ignore Point{D,T} = FrechetDist.cg.point.Point{D,T}
@ignore Dist =  FrechetDist.cg.point.Dist
@ignore Points =  FrechetDist.cg.polygon.Points

######################################################################
# Misc
######################################################################

@inline function  gid(P::Point{D,T},
                       r::T)::Point{D,Int} where{D,T}
    x = MVector{D,Int}(undef);

    for i ∈ 1:D
        x[i] = floor( Int, P[ i ] / r );
    end

    return  Point{D,Int}( x );
end


######################################################################
# Weighted point set
######################################################################
mutable struct WPoints{D,T}
    orig_PS::Vector{Point{D,T}};    # Original point setup
    PS::Vector{Int};                # The point set -- locations are
                                    # pointers to original points.
    W::Vector{Int};                 # Weights of points
end

@inline function  Base.length( P::WPoints{D,T} ) where {D,T}
    return  length( P.PS );
end


# An empty weighted point set, using the same ground set for defining
# the point set
function  WPoints_empty( G::WPoints{D,T} )  where{D,T}
    return  WPoints( G.orig_PS, Vector{Int}(), Vector{Int}() );
end

function  WPoints( _PS::Vector{Point{D,T}} )  where{D,T}
    PS = [i for i ∈ 1:length(_PS) ]
    W = fill(1, length( PS ) );

    return  WPoints( _PS, PS, W );
end

@inline function  Base.getindex(P::WPoints{D,T}, i::Int) where {D,T}
    return   P.orig_PS[ P.PS[ i ] ];
end


@inline function  weight(P::WPoints{D,T}, i::Int) where {D,T}
    return   P.W[ i ];
end


function   add_weight_to_point( P::WPoints{D,T}, i::Int, w::Int ) where {D,T}
    P.W[ i ] += w;
end

function  oindex( P::WPoints{D,T}, i ) where {D,T}
    return  P.PS[ i ];
end

function  add_point!( P::WPoints{D,T}, oind::Int, w::Int ) where {D,T}
    push!( P.PS, oind );
    push!( P.W , w );
    @assert( length( P.PS ) == length( P.W ) );
    return  length( P.PS );
end


######################################################################
# Linear time algorithm using hashing - using vectors to implement
# open hashing - probably a terrible idea.
######################################################################


mutable struct WGrid{D,T}
    P::WPoints{D,T}
    ℓ::Float64  # Side length
    cells::Dict{Point{D,Int},Vector{Int}}
end


"""
    store_value!

    An insertion function for dictionary where a value is an array of values.
"""
@inline function store_value!(dict::Dict{Point{D,Int},Vector{Int}},
                            key::SVector{D,Int}, value::Int) where {D}
    if haskey(dict, key)
        push!(dict[key], value)
    else
        dict[key] = [value]
    end
end

function   WGrid_init( Ground::WPoints{D,T}, ℓ::T )::WGrid{D,T}  where{D,T}
    dict = Dict{Point{D,Int},Vector{Int}}();

    return  WGrid( WPoints_empty( Ground ), ℓ, dict );
end

function   packing_add_point( G::WGrid{D,T}, loc::Int, rad::T,
                              P::WPoints{D,T} ) where{D,T}
    p =  P[ loc ];
    trg = gid( p, G.ℓ )
    id_min = trg .- 1;
    id_max = trg .+ 1;

    for  cell ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        i = Point{D,Int}( Tuple( cell )... );
        if  ! haskey( G.cells, i )  continue  end
        list = G.cells[ i ];
        for  p_ind ∈ list
            if   Dist( G.P[ p_ind ], p ) < rad
                add_weight_to_point( G.P, p_ind, weight( P, loc ) );
                return
            end
        end
    end

    oloc = oindex( P, loc )
    i = add_point!( G.P, oloc, weight( P, loc ) );
    store_value!( G.cells, trg, i );
end


function  packing( P::WPoints{D,T}, rad::T ) where{D,T}
    G = WGrid_init( P, rad );
    for  i ∈ 1:length( P )
        packing_add_point( G, i, rad, P );
    end

    return  G.P;
end


function  packing( _P::Vector{Point{D,T}}, rad::T ) where{D,T}
    P = WPoints( _P );

    N = packing( P, rad );

    out = Vector{Int}();
    for  i ∈ 1:length(N)
        push!( out, oindex( N, i ) );
    end

    return  out;
end


############################################################################
############################################################################


function   grid_neighberhood( base::Point{D,Int}, dist::T, ℓ::T ) where {D,T}
    Δ = ceil(Int, dist / ℓ ) + 1;

    id_min = base .- Δ
    id_max = base .+ Δ

    #    ℓ = dist / sqrt( D )
    #n = length( P )
    #G = grid_store( P, ℓ, 1:n )
    nbr = Vector{Point{D,Int}}();
    for  _subc ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        subc = Point{D,Int}( Tuple( _subc )... )
        sum::Float64 = 0
        for i ∈ 1:D
            x::Float64 = max( abs( subc[ i ] - base[ i ] ) - 1, 0 );
            sum += x*x;
        end
        d = sqrt(sum) * ℓ;
        if  ( d > dist )
            continue;
        end
        push!( nbr, subc );
    end

    return  nbr;
end


@inline function  check_neighbors( G::WGrid{D,T}, cell, i_pnt,
                                   far, rad, nbr ) where{D,T}
    P = G.P;
    p = P[ i_pnt ];
    for  _subc ∈ nbr
        subc = cell + _subc;
        if  haskey( G.cells, subc )
            subl = G.cells[ subc ]
            for  j ∈ subl
                if  j == i_pnt
                    continue
                end
                if  Dist( p, P[ j ] ) <= rad
                    far[ i_pnt ] = far[ j ]= false
                    return
                end
            end
        end
    end
end

"""
    grid_store

Stor all the points in P[ rng ] in the grid G.
"""

function   grid_store( G, P::WPoints{D,T}, rng ) where{D,T}
    for i ∈ rng
        p = P[ i ]
        trg = gid( p, G.ℓ )
        add_point!( G.P, oindex( P, i ), P.W[ i ] )
        store_value!( G.cells, trg, length( G.P) );
    end
end

function  r_far( P::WPoints{D,T}, rad::T ) where{D,T}
    ℓ = rad / sqrt( D )
    n = length( P )
    G = WGrid_init( P, ℓ )

    nbr = grid_neighberhood( zero(Point{D,Int}), rad, ℓ );

    grid_store( G, P, 1:n )

    #Δ = ceil(Int, sqrt( D ) )

    far = fill( true, n )
    for  (cell, list) ∈ G.cells
        # If there are 2 or more points in the cell, then they are all near...
        if  length( list ) > 1
            for i ∈ list
                far[ i ] = false;
            end
            continue;
        end
        check_neighbors( G, cell, list[ 1 ], far,  rad, nbr );
    end

    return  far;
end

############################################################################
############################################################################


function  nn_dist( P, ind::Int )
    n = length( P );
    @assert( n > 1 );
    p = P[ ind ];
    min_ind = ( ind == 1 ) ? 2 : 1
    min_dist = Dist( p, P[ min_ind ] );
    for  j ∈ 1:n
        if  j == ind  continue  end

        d = Dist( p, P[ j ] )
        if  d < min_dist
            min_dist = d
            min_ind = j
        end
    end

    return  min_ind, min_dist
end

"""
    random_nn_dist

# Returned value

Returns a triple: i, j, dist.

i: the index of a random point in P,
j: The index of the closest point in P to P[i].
dist: Distance between P[i] and P[j].
"""
function  random_nn_dist( P )
    i = rand( 1:length( P ) );
    return  i, nn_dist( P, i )...;
end


function  max_weight( P::WPoints{D,T} ) where {D,T}
    i = argmax( P.W );
    return  i, P.W[ i ]...;
end



function  subset( P::WPoints{D,T}, keep ) where {D,T}
    @assert( length( P ) == length( keep ) );

    O = WPoints_empty( P );

    for  i ∈ 1:length( P.W )
        if   keep[ i ]
            add_point!( O, oindex( P, i ), P.W[ i ] )
        end
    end

    return  O;
end


function   min_ball( P::Vector{Point{D,T}}, i, k ) where {D,T}
    cen = P[ i ];

    arr = [Dist(cen, q) for q ∈ P ];
    d = partialsort!( arr, k )

    return  cen, d;
end


############################################################################
############################################################################
function   smallest_k_ball( _P::Vector{Point{D,T}}, k::Int ) where{D,T}
    P = WPoints( _P );

    while  true
        _, _, dist = random_nn_dist( P );

        N = packing( P, dist );
        mw = max_weight( N )[ 2 ]

        # If there is a packing point with weight ≥ k, then we can
        # throw away the far points
        if  ( mw >= k )
            far = r_far( P, dist );
            near = .!far;
            if ( sum( near ) < length( P ) )
                P = subset( P, near );
                continue;
            end
        end

        N_4 = packing( P, 4*dist );
        ind, mw_4 = max_weight( N_4 )
        if  ( mw_4 >= k )
            return  min_ball( _P, oindex( P, ind ), k );
        end

        P = N;
    end
end


############################################################################
############################################################################

function  draw_two_sets( P, N, filename )
    # Output to file...
    c,cr,_ = cairo_setup( filename, [P], true );

    n = length( P );
    set_source_rgb(cr, 0.8, 0.2, 0.2) # An nice red color
    draw_points( cr, Points( N ), 1.0/(sqrt(n)) );
    Cairo.stroke( cr );

    set_source_rgb(cr, 1.0, 0.0, 0.8);
    draw_points( cr, Points( P ), 1/(16*sqrt(n)) );
    Cairo.stroke( cr );

    finish( c );
    println( "     Created: ", filename );
end


function  test_far_points( P, rad )
    far = r_far(WPoints( Points( P ) ), rad );

    N = Polygon2F();
    for  i ∈ eachindex( P )
        if  ( far[ i ] )
             push!( N, P[ i ] );
        end
    end

    draw_two_sets( P, N, "out/far_points.pdf" );
end


function  test_packing( P, rad )
    net = packing( Points( P ), rad );

    N = Polygon2F();
    for  i ∈ net
        push!( N, P[ i ] );
    end

    draw_two_sets( P, N, "out/net.pdf" );
end



function     test_smallest_ball( P, k )
    cen, rad = smallest_k_ball( Points( P ), k );

    c,cr,_ = cairo_setup( "out/smallest_disk.pdf", [P], true );

    set_source_rgb(cr, 0.8, 1.0, 0.2) # An nice red color
    Cairo.arc(cr, cen[1], cen[2], rad, 0, 2*pi)
    Cairo.fill(cr)
    Cairo.stroke( cr );

    n = length( P );
    set_source_rgb(cr, 0.8, 0.2, 0.2) # An nice red color
    draw_points( cr, Points( P ), 0.1/(sqrt(n)) );
    Cairo.stroke( cr );
    finish( c );
end

function (@main)(ARGS)
    if  length( ARGS ) != 2
        println( "Usage:\n\t"
                 * "packing.jl [n] [rad]\n\n" );
        return  -1;
    end

    D = 2;
    n =  parse(Int64, replace( ARGS[ 1  ], "," => "" ) );
    rad =  parse(Float64, replace( ARGS[ 2  ], "," => "" ) );

    println( "n: ", n );
    println( "rad: ", rad );

    println( "\n\n--------------------\n\n" );
    println( "Generating input & Copying..." );
    P = Polygon_random_gaussian( D,Float64, n );

    test_far_points( P, rad );
    test_packing( P, rad );

    test_smallest_ball( P, 10 )



    return  0;
end

#! /usr/bin/env  julia

# shortcut.jl:
#
# Approximate the optimal shortcut for a give data file.
#

push!(LOAD_PATH, pwd()*"/src/")
push!(LOAD_PATH, pwd()*"/src/cg/")

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using VirtArray;
using Printf;
using DataFrames
using PrettyTables

include( "graphics.jl" )

import BBT
import WSPD
using Cairo, Colors


# Proximity graph is a directed graph. The nodes are points, edges are
# between points. Each edge has quite a bit of

struct GGEdge
    src::Int
    dst::Int

    dist::Float64;
    diam::Float64; # The size of the "target" region the edge leads to.
end

mutable  struct GeometricGraph{D,T}
    PS::Polygon{D,T}         # The points
    E::Vector{GGEdge}    # The edges
    edges_start::Vector{Int}    # The start index of the outgoing edges of a vertex
end

import Base: isless

function  isless(a::GGEdge, b::GGEdge)
    if  a.src < b.src
        return  true
    end
    if  a.src > b.src
        return  false
    end

    return  a.dist < b.dist
end

#function sort_edges!( Edges::Vector{GGEdge} )
#end

function  Geometric_Graph( P::Polygon{D,T}, edges::Vector{GGEdge} ) where{D,T}
    sort!( edges );
    ES = zeros( length( P ) );
    for  i ∈ 1:length(edges )
        e = edges[ i ];
        #println( e.src );
        if  ES[ e.src ] == 0
            ES[ e.src ] = i;
        end
    end

    return  GeometricGraph{D,T}( P, edges, ES );
end


function  GeometricGraph_via_WSPD( P::Polygon{D,T}, ε::Float64 ) where {D,T}
    W = WSPD.init( P,  ε / 2  );
    WSPD.expand( W );
    tree = W.tree;
    edges = Vector{GGEdge}();

    function  add_edges( src::Int, node )
        for  i ∈ node.r
            i_p = orig_index( tree.PS, i );

            # Create the edge P[l_i] -> P[ i_p]
            e = GGEdge( src, i_p, Dist( P[ src ], P[ i_p ] ), node.diam );
            push!( edges, e );
        end
    end

    for  pair ∈ W.finals
        i_l,i_r = WSPD.reps_orig_indexes( W, pair );

        right = pair.right;
        add_edges( i_l, pair.right );
        add_edges( i_r, pair.left );
    end
    G = Geometric_Graph( P, edges );

    return   G;
end


function (@main)(ARGS)
    if  length( ARGS ) != 1
        println( "Usage:\n\t"
                 * "nn_graph.jl [n]\n\n" );
        return  -1;
    end

    D = 2;
    n =  parse(Int64, ARGS[ 1  ] );

    P = Polygon_random( D,Float64, n );
    ε = 0.5;

    G = GeometricGraph_via_WSPD( P, ε );

    println( "# edges: ", length( G.E ) );
    return  0;
end


####################################################################3

#export  BBTree_build;

#end

module  WSPD


push!(LOAD_PATH, pwd()*"/cg/")
push!(LOAD_PATH, pwd() )

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using FrechetDist.cg.bbox;
using VirtArray;
using Cairo, Colors

#include("./VirtArray.jl")
#include("./BBT.jl")

import BBT
using DataStructures

mutable struct WSPair{D,T}
    left::BBT.Node
    right::BBT.Node

    dist::T # Distance between the boxes.
    max_sides_diam::T
    diam_ub::T  # Upper bound on the maximum diameter of any pair of
                # points in this pair.
end

function  WPDPair_init( left::BBT.Node{D,T}, right::BBT.Node{D,T} )  where {D,T}
    @assert( !isnothing( left ) );
    @assert( !isnothing( right ) );

    dist = bbox.dist( left.bb, right.bb );
    max_sides_diam = max( bbox.diam( left.bb ), bbox.diam( right.bb ) );
    diam_ub = bbox.max_dist( left.bb, right.bb );

    return  WSPair{D,T}( left, right, dist, max_sides_diam, diam_ub );
end

const MAX_SEP = 100000000.0
function  separation( pair::WSPair{D,T} ) where {D,T}
    if  pair.dist <= 0.0
        #if  pair.max_sides_diam <= 0.0
        #    return  0.0;
        #end
        return  MAX_SEP   # Or should I return +∞
    end
    return  pair.max_sides_diam / pair.dist;
end


#mutable struct  WSPDOrder{T} <: Base.Order.Ordering
#
#end

mutable struct  WSPDOrder{D,T} <: Base.Order.Ordering
    dummy::Int64
end

#function  lt(o::EventsOrder, a::Int64, b::Int64 )
#    return  isless( o.values[ a ], o.values[ b ] );
#end

import Base.Order.lt
function  lt( o::WSPDOrder{D,T}, a::WSPair{D,T}, b::WSPair{D,T} ) where {D,T}
    o.dummy += 1;
    return  a.diam_ub > b.diam_ub;
end

PairInt = Tuple{Int,Int};
mutable struct PD{D,T}
    pnts::Polygon{D,T};
    finals::Vector{WSPair{D,T}};
    curr_active::Vector{WSPair{D,T}};
    heap::BinaryHeap{WSPair{D,T}};
    tree::BBT.Tree{D,T}
    sep::T;  # Desired quality of separation
    hash_pairs::Dict{PairInt, Bool};
end


@inline function  get_id( u::BBT.Node, v::BBT.Node )::PairInt
    x, y = u.id, v.id;

    x = min( u.id, v.id );
    y = max( u.id, v.id );

    return  (x,y)
end


@inline function  has_pair( W::PD{D,T}, pair::PairInt ) where {D,T}
    return  haskey( W.hash_pairs, pair );
end


function  push_pair( W::PD{D,T}, u::BBT.Node{D,T}, v::BBT.Node{D,T}
                   ) where {D,T}
    id = get_id( u, v );
    if  has_pair( W, id )
        #println( "REJECTED:", id );
        return
    end;

    #println( id );
    p_a = WPDPair_init( u, v );
    push!( W.heap, p_a );
    W.hash_pairs[ id ] = true;

    return  p_a;
end

"""
    top_refine

    Takes the top of the active pairs, and refine it.
"""
function  top_refine( W::PD{D,T} ) where {D,T}
    if ( isempty( W.heap ) ) return  end

    top = pop!( W.heap );
    BBT.Tree_refine_node( W.tree, top.left );
    BBT.Tree_refine_node( W.tree, top.right );

    l_diam = bbox.diam( top.left.bb );
    r_diam = bbox.diam( top.right.bb );

    # Is this pair already between points?
    if  ( ( r_diam == 0.0 )  &&  ( l_diam == 0.0 ) )
        # A silly self pair of a point with itself.
        if  top.dist <= 0.0
            return;
        end

        push!( W.finals, top );
        return;
    end

    f_refine_left = (l_diam >= r_diam );

    if  f_refine_left
        push_pair( W, top.left.left, top.right );
        push_pair( W, top.left.right, top.right );
    else
        push_pair( W, top.left, top.right.left );
        push_pair( W, top.left, top.right.right );
    end
end


"""
    top_delete

    Takes the top of the active pairs, and deletes it.
"""
function  top_delete( W::PD{D,T} )  where{D,T}
    if ( isempty( active ) ) return  end

    top = pop!( W.heap );
end


"""
    top_finalize

Takes the top of the active pairs, and move it to the generated list of pairs.
"""
function  top_finalize( W::PD{D,T} )  where  {D,T}
    if ( isempty( W.heap ) ) return  end
    top = pop!( W.heap );
    push!( W.finals, top );
end


function  init( _pnts::Polygon{D,T}, _sep::T ) where {D,T}
    tree = BBT.Tree_init( _pnts );

    PairT = WSPair{D,T}
    #order = Base.Order.ForwardOrdering{PairT}( compare_pairs );
    order =  WSPDOrder{D,T}( 0 );

    curr_active = Vector{PairT}();
    heap = BinaryHeap{PairT}( order, curr_active );
    W = PD( _pnts, Vector{WSPair{D,T}}(),
        curr_active, heap, tree, _sep,
        Dict{PairInt, Bool}() );

    pair = WSPair{D,T}( tree.root, tree.root, 0.0, bbox.diam( tree.root.bb ),
                        bbox.diam( tree.root.bb ) );
    push!( W.heap, pair );

    return  W;
end


"""
    expand( W )

Construct the whole WSPD.

"""
function   expand( W::PD{D,T} ) where {D,T}
    while  ( !isempty( W.heap ) )
        top = first( W.heap );
        #println( "    Δ: ", top.max_diam, "    [", separation( top ), "]:",
        #    length( W.finals ) );
        if  ( separation( top ) > W.sep )
            top_refine( W );
            continue;
        end
        #println( top.left.r, " × ", top.right.r );
        top_finalize( W );
    end
end

"""
     expand_bichromatic

expand pairs, but only if they are not fully contained either in
1:index_cut range or the index_cut+1:end range.  in the original
array. Gives an easy way to get WSPD for bichromatic pairs.

"""
function   expand_bichromatic( W::PD{D,T}, index_cut ) where {D,T}
    while  ( !isempty( W.heap ) )
        top = first( W.heap );

        i_min = min_orig_index( W, top );
        i_max = max_orig_index( W, top );
        f_boring = ( ( ( i_min <= index_cut )  &&  (i_max <= index_cut ) )
                     || ( ( i_min > index_cut )  &&  (i_max > index_cut ) ) );
        if   f_boring
            pop!( W.heap );
            continue;
        end

        if  ( separation( top ) > W.sep )
            top_refine( W );
            continue;
        end
        #println( top.left.r, " × ", top.right.r );
        top_finalize( W );
    end
end

function   get_top( W::PD{D,T} ) where {D,T}
    @assert( ! isempty( W.heap ) );

    return  first( W.heap );
end

function   top_diam_ub( W::PD{D,T} ) where {D,T}
    @assert( ! isempty( W.heap ) );
    if isempty( W.heap )
        return  zero( T );
    end

    return  first( W.heap ).diam_ub;
end

function  get_reps( W::PD{D,T}, pair::WSPair{D,T} ) where {D,T}
    @assert( ! isempty( pair.left.r ) )
    @assert( ! isempty( pair.right.r ) )
    return  W.tree.PS[ first( pair.left.r ) ], W.tree.PS[ first( pair.right.r ) ]
end

function  reps_orig_indexes( W::PD{D,T}, pair::WSPair{D,T} ) where {D,T}
    @assert( ! isempty( pair.left.r ) )
    @assert( ! isempty( pair.right.r ) )

    l_i = first( pair.left.r )
    r_i = first( pair.right.r )

    #println( typeof( W.tree.PS ) );

    #println( orig_index( W.tree.PS, l_i ) );
    return  orig_index( W.tree.PS, l_i ),orig_index( W.tree.PS, r_i )
end


function  min_orig_index( W::PD{D,T}, pair::WSPair{D,T} ) where {D,T}
    @assert( ! isempty( pair.left.r ) )
    @assert( ! isempty( pair.right.r ) )

    return  min( BBT.get_min_orig_index( W.tree, pair.left ),
                 BBT.get_min_orig_index( W.tree, pair.right ) );
end

function  max_orig_index( W::PD{D,T}, pair::WSPair{D,T} ) where {D,T}
    @assert( ! isempty( pair.left.r ) )
    @assert( ! isempty( pair.right.r ) )

    return  max( BBT.get_max_orig_index( W.tree, pair.left ),
                 BBT.get_max_orig_index( W.tree, pair.right ) );
end

function  get_orig_ranges( W::PD{D,T}, pair::WSPair{D,T} ) where {D,T}
    @assert( ! isempty( pair.left.r ) )
    @assert( ! isempty( pair.right.r ) )

    l_min = BBT.get_min_orig_index( W.tree, pair.left )
    l_max = BBT.get_max_orig_index( W.tree, pair.left )
    r_min = BBT.get_min_orig_index( W.tree, pair.right )
    r_max = BBT.get_max_orig_index( W.tree, pair.right )
    
    return  l_min:l_max,r_min:r_max
end


### The exports...

export WPD
export WPDPair

export top_delete, top_finalize, top_delete
export expand, init;

export expand_bichromatic

export top_diam_ub
export get_top
export top_refine

export get_reps
export reps_orig_indexes

export min_orig_index
export max_orig_index
export get_orig_ranges

end

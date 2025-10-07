####################################################
# BBT - BoundingBoxTree
#
# A tree that stores points in a hierarchical axis aligned bounding boxes.
####################################################

module  BBT

#push!(LOAD_PATH, pwd()*"/cg/")
#push!(LOAD_PATH, pwd() )

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using FrechetDist.cg.bbox;

#include("./VirtArray.jl")
#include( "VirtArray.jl" );
using VirtArray;
using Cairo, Colors

#using Base: reindex
IntRange = UnitRange{Int};

mutable struct Node{D,T}
    bb::BBox{D,T}
    r::IntRange
    left::Union{Nothing, Node}
    right::Union{Nothing, Node}
    f_leaf::Bool;
    diam::T;
    id::Int
    min_index::Union{Nothing, Int}
    max_index::Union{Nothing, Int}
end

mutable struct Tree{D,T}
    PS::VArray{Point{D,T}};
    root::Union{Nothing, Node}
    id_counter::Int
end

########################################################################


function  get_min_max_orig_index( tree::Tree, node::Node,
    min_max_func, field_name::Symbol )::Int

    if  ! isnothing( getfield( node, field_name ) )
        return  getfield( node, field_name );
    end
    if   node.f_leaf
        setfield!( node, field_name,  orig_index( tree.PS, first( node.r ) ) );
        return  getfield( node, field_name );
    end

    # We could potentially expand the node, and then use the
    # following, but this seems somewhat wasteful. For the time being,
    # if we reached a current un-expanded node, we just compute the
    # minimum using brute force.

    # node_expand( node, tree );

    if  ( ( ! isnothing( node.left ) )  &&  ( ! isnothing( node.right ) ) )
        v_l = get_min_max_orig_index( tree, node.left, min_max_func, field_name );
        v_r = get_min_max_orig_index( tree, node.right, min_max_func, field_name ) 
        val = min_max_func( [ v_l v_r ] );
    else
        val = min_max_func( [orig_index( tree.PS, i ) for i ∈ node.r] )
    end
        
    setfield!( node, field_name, val );
    return  getfield( node, field_name );
end


function  get_min_orig_index( tree::Tree, node::Node )::Int
    return  get_min_max_orig_index( tree, node, minimum, :min_index );
end

function  get_max_orig_index( tree::Tree, node::Node )::Int
    return  get_min_max_orig_index( tree, node, maximum, :max_index );
end

    #=
function  Node_get_min_index_old( tree::Tree, node::Node )::Int
    if  ! isnothing( node.min_index )
        return  node.min_index
    end
    if   node.f_leaf
        node.min_index = orig_index( tree.P, first( r ) );
        return  node.min_index;
    end

    # We could potentially expand the node, and then use the
    # following, but this seems somewhat wasteful. For the time being,
    # if we reached a current un-expanded node, we just compute the
    # minimum using brute force.

    # node_expand( node, tree );

    if  ( ( ! isnothing( node.left ) )  &&  ( ! isnothing( node.right ) ) )
        node.min_index = minimum( [Node_get_min_index( tree, node.left ),
                                   Node_get_min_index( tree, node.right ) ] );
        return  node.min_index;
    end

    node.min_index = minimum( [orig_index( tree.P, i ) for i ∈ r] );

    return  node.min_index;
end

=#




function  node_init( tree::Tree{D,T}, range::IntRange )  where {D,T}
    bb = BBox{D,T}();
    tree.id_counter += 1;

    #println( "NID: ", tree.id_counter );
    bbox.bound( bb, tree.PS[ range ] );
    return  Node( bb, range, nothing, nothing, false, bbox.diam( bb ),
                     tree.id_counter, nothing, nothing );
end

function  Tree_init_inner( p::Polygon{D,T} )::Tree{D,T} where {D,T}
    varr = VArray( Points( p ) );
    tree = Tree( varr, nothing, 0 );

    tree.root = node_init( tree, 1:length(p) );

    return  tree;
end

function  Tree_init_inner( _p::Vector{Point{D,T}} )::Tree{D,T} where {D,T}
    p = Polygon{D,T}( _p );
    return  Tree_init_inner( p );
end



function  is_identical( arr::VArray{T}, r )::Bool where {T}
    if  ( length( r ) <= 1 )   return  true;   end

    s = first( r )
    for  i  in  s+1:last(r)
        if  arr[ s ] != arr[ i ]
            return  false;
        end
    end

    return  true;
end


function pnt_varray_partition!( P::VArray{Point{D,T}}, r::IntRange,
                                dim, cutoff )  where {D,T}
    @assert( ( 1 <= first( r ) )  &&  ( last(r) <= length( P ) ) );

    len = top = last( r );
    @inbounds for  i  in  r
        (i > top)  &&  break;
        while  ( P[ i ][ dim ] > cutoff )
            swap!( P, i, top );
            #P[ i ], P[ top ] = P[ top ], P[ i ];
            top = top - 1;
            (i > top)  &&  break;
        end
    end
    return  first(r):top, (top + 1):len;
end


function  node_split( node::Node{D,T}, tree::Tree{D,T} )  where {D,T}
    if  ( ( ! isnothing( node.left ) )  &&  ( ! isnothing( node.right ) ) )
        #println( "Already split?" );
        return;
    end

    # If a leaf there is nothing to split.
    if  node.f_leaf  return  end


    #println( "node_split" );
    # covers the case range isa 1
    if  is_identical( tree.PS, node.r )
        #println( "Leaf???" );
        node.f_leaf = true;
        return;
    end

    #println( "Doing split!" );
    w, dim = findmax( i->cg.width( node.bb, i ), 1:D );

    cutoff = bbox.middle( node.bb, dim );
    r_l, r_r = pnt_varray_partition!( tree.PS, node.r, dim, cutoff )
    @assert( ( length( r_l ) > 0 )  &&  ( length( r_r ) > 0 ) );

    node.left = node_init( tree, r_l );
    node.right = node_init( tree, r_r );
    #println( "L NID: ", node.left.id )
    #println( "R NID: ", node.right.id )

    return  node;
end


function  node_expand( v::Node{D,T}, tree::Tree{D,T} ) where {D,T}
    v.f_leaf  &&  return;

    node_split( v, tree );

    #    v.f_leaf  &&  return;
    #
#    ( v.left != nothing )   &&  node_expand( v.left, tree );
#    ( v.right != nothing )  &&  node_expand( v.right, tree );
end

function Tree_init( PS::Polygon{D,T} )::Tree{D,T} where {D,T}
    tree = Tree_init_inner( PS );

    #println( typeof( tree ) );
    #println( "tree created" );

    node_expand( tree.root, tree );

    return  tree;
end

function  fully_expand( node::Node{D,T}, tree::Tree{D,T} ) where {D,T}
    if   node.f_leaf  return  end;

    if ( isnothing( node.left )  ||  isnothing( node.right ) )
        #println( "   mode expand? " );
        node_expand( node, tree );
    end
    if  ( ! isnothing( node.left ) )
        fully_expand( node.left, tree );
    end
    if  ( ! isnothing( node.right ) )
        fully_expand( node.right, tree );
    end
end

function  Tree_fully_expand( tree::Tree{D,T} ) where {D,T}
    fully_expand( tree.root, tree );
end

function   bbox_draw( context, bb, color )
    set_source(context, color)

    bl = bbox.bottom_left( bb );
    w = cg.width( bb );
    h = cg.height( bb );

    rectangle(context, bl[1], bl[2], w, h );
    fill( context );
end


function   node_draw( context, node, level, range )
    if  node.f_leaf  return  end;

    if  level ∈ range
        # 0.1 alpha means 10% opaque, 90% transparent
        yellow_transparent = coloralpha(colorant"yellow", 0.1)

        set_source(context, yellow_transparent)

        bb = bbox.expand( node.bb, 1.01 );

        bl = bbox.bottom_left( bb );
        w = width( bb );
        h = height( bb );

        # Draw the rectangle
        rectangle(context, bl[1], bl[2], w, h );
        fill_preserve(context)

        # Draw the blue frame
        set_source(context, colorant"blue")
        set_line_width(context, 1.0) # Set the line width for the frame
        stroke(context)
    end
    if  ( ! isnothing( node.left ) )
        node_draw( context, node.left, level+1, range );
    end
    if  ( ! isnothing( node.right ) )
        node_draw( context, node.right, level+1, range );
    end
end

function  depth( node::Node{D,T} )::Int where  {D,T}
    if  ( isnothing( node ) )  return  0  end;
    if  node.f_leaf  return 1 end;

    return 1 + max(  depth( node.left ), depth( node.right ) )
end

function Tree_draw( tree::Tree{D,T}, filename::String ) where{D,T}
    bb = bbox.expand( tree.root.bb, 1.3 );
    #bbox.print( bb );

    surface = CairoPDFSurface(filename, bbox.max( bb, 1 ),
        bbox.max( bb, 2 ) )
    context = CairoContext(surface)

    bbox_draw( context, bb, colorant"lightblue" );

    d = depth( tree.root );
    node_draw( context, tree.root, 0, 0:d );

    # Set the fill color to yellow with 90% transparency

    show_page(context)

    # Now draw levels of the tree...
    for  i ∈ (d-1):-1:0
        bbox_draw( context, bb, colorant"lightblue" );
        node_draw( context, tree.root, 0, i:i+1 );
        show_page( context )
    end

    finish(surface)
    println( "\n" * "Created file: $filename")
end


function  Tree_refine_node( tree, node )
    node_split( node, tree );
end

# Export list...

export  Node,  Tree;
export  Tree_draw, Tree_init, Tree_fully_expand
export  Tree_refine_node

export   get_min_orig_index;
export   get_max_orig_index;

end  # Module
####################################################################3

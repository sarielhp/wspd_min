#! /usr/bin/env  julia

# shortcut.jl:
#
# Approximate the optimal shortcut for a give data file.
#

push!(LOAD_PATH, pwd()*"/src/")
push!(LOAD_PATH, pwd()*"/src/cg/")

using Statistics
using FrechetDist
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

Point3F = Point{3,Float64};
DiskF = Point3F;
        

function (@main)(ARGS)
    if  length( ARGS ) != 1
        println( "Usage:\n\t"
                 * "nn_graph.jl [n]\n\n" );
        return  -1;
    end

    D = 2;
    n =  parse(Int64, ARGS[ 1  ] )

    Q = Polygon_random( 2, Float64, n )
    ε = 0.5;

    v = Statistics.mean( Points( Q ) )
    
    P = Polygon2F( filter( p -> Dist(p, v) > 0.1, Points( Q ) ) )
    
    #P = Q
    println( typeof( P ) );
    N = sort( Points( P ), by = p -> Dist(p, v ) )
    push!( P, Point2F( v ) )

    disks = Vector{DiskF}();
    reps = Polygon2F();
    cens = Polygon2F();
    α = 4;
    while length( N ) > 0
        p = popfirst!( N )

        push!( reps, p );
        τ = 1/(α^2 -1);
        cen = p + τ * ( p - v ); 
        #println( "p: ", p, "  v: ", v, "  cen: ", cen );
        rad = α / (α^2 - 1) * Dist( p, v );
        
        d = DiskF( cen[ 1 ], cen[ 2 ], rad )
        #println( d );

        push!( cens, cen );
        push!( disks, d );
        filter!( p -> Dist( cen, p ) > rad, N )
        for z in N
            if Dist( cen, z ) <= rad
                println( "ZZZZ :", z )
            end
        end
    end

    c,cr,bb = cairo_setup( "out.pdf", [P], true )

    #draw_points( cr, Points( P ) )
    
    #Cairo.show_page( cr );

    m_p = Polygon2F();
    push!( m_p, Point2F(v) );


    set_source_rgb( cr, get_color_rgb( 1 )... );
    draw_circles( cr, disks )
    set_source_rgb( cr, get_color_rgb( 5 )... );
    draw_points( cr, Points( m_p ), 0.01 )
    set_source_rgb( cr, get_color_rgb( 2 )... );
    draw_points( cr, Points( cens ), 0.002 )
    set_source_rgb( cr, get_color_rgb( 4 )... );
    draw_points( cr, Points( reps ), 0.003 )
    #set_source_rgb( cr, get_color_rgb( 5 )... );
    #draw_points( cr, Points( reps ), 0.002 )

    Cairo.show_page( cr );
    set_source_rgb( cr, get_color_rgb( 5 )... );
    draw_points( cr, Points( m_p ), 0.01 )
    set_source_rgb( cr, get_color_rgb( 4 )... );
    draw_points( cr, Points( reps ), 0.003 )

    println( "# points in P: ", length( P ) );

    
    #=Cairo.show_page( cr );
    set_source_rgb( cr, get_color_rgb( 1 )... );
    draw_points( cr, Points( P ) )
    =#
    #=
    Cairo.show_page( cr );

    draw_points( cr, Points( P ) )
    set_source_rgb( cr, get_color_rgb( 2 )... );
    draw_points( cr, Points( reps ) )
    set_source_rgb( cr, get_color_rgb( 3 )... );
    draw_circles( cr, disks )
    =#
    

    finish( c )

    return  0;
end


####################################################################3

#export  BBTree_build;

#end

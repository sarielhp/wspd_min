#! /usr/bin/env  julia

# diameter_smf.jl:
#
# Compute the diameter of an smf data file.
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

using BBT
using WSPD
using Cairo, Colors


function  rt_str( rt )
    return  @sprintf( "%.6f", rt );
end

function  rt_str_2( rt )
    return  @sprintf( "%.2f", rt );
end

function  approx_diameter( P::Polygon{D,T}, ε::Float64  ) where {D,T}
    if  length( P ) <= 1
        return  0.0
    end
    c = 1.0 + ε;
    W = WSPD_init( P,  ε / 2  );
    curr = Dist( P[ 1 ], P[ 2 ] );
    while  ( WSPD_top_diam_ub( W ) > c * curr )
        top = WSPD_get_top( W );
        p,q = WSPD_get_reps( W, top );
        curr = max( curr, Dist( p, q ) );
        WSPD_top_refine( W );
    end

    return  curr
end

function (@main)(ARGS)
    if  length( ARGS ) != 2
        println( "Usage:\n\t"
                 * "diameter_smf.jl [eps] [smf.file]\n\n" );
        return  -1;
    end
#    println( ARGS[ 1 ] );
    eps =  parse(Float64, ARGS[ 1  ] );
    P = read_smf_file( ARGS[ 2 ] );

    approx_diameter( P, 2.0 );
    t = @timed diam = approx_diameter( P, eps );
    println( t.time );

    @printf( "%.6f-approx diameter: %15.6f   Time: %11.6f    N: %7d\n",
             eps, diam, t.time, length( P ) );
#    println( eps"-approx diameter: ", diam,  " time: " );
    return  0;
end


####################################################################3

#export  BBTree_build;

#end

#! /usr/bin/env  julia
#
# sphere_diameter.jl
#
# Compute the diameter of a random point set on the sphere of various
# dimensions.
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

function  BBT_test( P, N )
    tree = BBTree_init( P );
    #println( "Tree initialized..." );
    BBTree_fully_expand( tree );
    BBTree_draw( tree, "test.pdf" );
end

function  exact_diameter( P::Polygon{D,T}  ) where {D,T}
    curr::Float64 = 0;

    for  i ∈ 1:length(P)-1
        for j ∈ i+1:length(P)
            curr = max( curr, Dist( P[i], P[j] ) );
        end
    end
    return  curr;
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

function  rt_str( rt )
    return  @sprintf( "%.6f", rt );
end

function  rt_str_2( rt )
    return  @sprintf( "%.2f", rt );
end

function  test_diameter( P, N::Int64, D )
    mult!( P, 900.0 );
    shift!( P, Point{D,Float64}( 200.0, 200.0, zeros( D - 2 )... ));

    t_approx = @timed approx_diam = approx_diameter( P, 0.1 );

    t_exact = @timed exact_diam = exact_diameter( P );

    @printf( "D: %d  N: %8d  diam ≈ %10g  Exact: %10g    T≈ %10.4f  ET= %10.4f\n",
        D, N, approx_diam, exact_diam, t_approx.time, t_exact.time );
    new_row_df = DataFrame(
        dimension=D,
        N=N,
        approx_runtime=rt_str(t_approx.time),
        exact_runtime=rt_str( t_exact.time ),
        approx_ratio=rt_str_2( approx_diam/exact_diam )
#        diam=exact_diam,
    )

    return  new_row_df
end


function  diam_test_sphere( D, iters, filename )
    # Force recompilation
    begin
        P = Polygon_random_sphere( D, Float64, 40 );
        approx_diameter( P, 0.1 );
        exact_diameter( P );
    end

    df = DataFrame();
    for  i ∈ 1:iters
        N = 2^i
        P = Polygon_random_sphere( D, Float64, N );
        new_row = test_diameter( P, N, D );
        #println( new_row );
        append!(df, new_row )
        #println( df );
    end

    open(filename, "w") do fl
        pretty_table( fl, df, header=[ "Dim", "N","RT Approx","RT Exact", "Approx"],
            alignment=:r,
            backend = Val(:markdown)  );
    end
    println( "Created: ", filename );
end

function (@main)(ARGS)
    mkpath( "results" );

    diam_test_sphere( 2, 19, "results/2d_sphere.md" )
    diam_test_sphere( 3, 19, "results/3d_sphere.md" )
    diam_test_sphere( 4, 19, "results/4d_sphere.md" )
    diam_test_sphere( 5, 19, "results/5d_sphere.md" )
    diam_test_sphere( 6, 19, "results/6d_sphere.md" )
    return  0;
end


####################################################################3

#export  BBTree_build;

#end

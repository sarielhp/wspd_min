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


function  rt_str( rt )
    return  @sprintf( "%.6f", rt );
end

function  rt_str_2( rt )
    return  @sprintf( "%.2f", rt );
end

function  approx_shortcut( P::Polygon{D,T}, ε::Float64  ) where {D,T}
    @assert( length( P ) > 1 );

    if  length( P ) <= 1
        return  0.0
    end

    plens = Polygon_prefix_lengths( P );

    c = 1.0 + ε;
    W = WSPD.init( P,  ε / 2  );
    WSPD.expand( W );

    max_quality = -1.0;
    p_i = q_i = 1;
    for  pair ∈ W.finals
        #i_min = WSPD.min_orig_index( W, pair );
        #i_max = WSPD.max_orig_index( W, pair );
        l_i,r_i = WSPD.reps_orig_indexes( W, pair );
        #println( "M: ", i_min, "/", i_max, "   ", l_i, "/", r_i );
        l = abs( plens[ l_i ] - plens[ r_i ] )
        d = Dist( P[ l_i ], P[ r_i ] );

        if  d == 0.0   continue  end

        quality = l / d
        if quality > max_quality
            p_i = min( l_i, r_i);
            q_i = max( l_i, r_i );
            max_quality = quality
        end
    end

    return  p_i, q_i, max_quality;
end


function  approx_bichromatic_shortcut( P::Polygon{D,T},
                                       prefix::Int,
                                       ε::Float64, min_quality::Float64 ) where {D,T}
    @assert( length( P ) > 1 );

    if  length( P ) <= 1
        return  1,1,1.0
    end

    plens = Polygon_prefix_lengths( P );

    c = 1.0 + ε;
    W = WSPD.init( P,  ε / 2  );
    WSPD.expand_bichromatic( W, prefix );

    max_quality = -1.0;
    p_i = q_i = length( P ) + 1;
    for  pair ∈ W.finals
        i_min = WSPD.min_orig_index( W, pair );
        i_max = WSPD.max_orig_index( W, pair );
        if   ( ( i_min > prefix )  ||  ( i_max <= prefix ) )
            continue
        end
        #l_i,r_i = WSPD.reps_orig_indexes( W, pair );
        #println( "M: ", i_min, "/", i_max, "   ", l_i, "/", r_i );
        l = abs( plens[ i_max ] - plens[ i_min ] )
        d = Dist( P[ i_max ], P[ i_min ] );

        if  d == 0.0   continue  end

        quality = l / d
        if quality > min_quality
            if  ( i_min < p_i )
                p_i = i_min;
                q_i = i_max;
                max_quality = quality;
            end
        end
    end

    return  p_i, q_i, max_quality;
end


function repeated_shortcut( P );

    P = Polygon_sample_uniformly( P, 1000 );
    println( "N = ", length( P ) );
    approx_shortcut( P, 1.0 );
    list = VecPolygon2F();
    for  i ∈ 1:40
        push!( list, P );
        t = @timed  v_i,v_j,quality = approx_shortcut( P, eps );
        println( i, " : ", quality );
        if  ( quality < 1.250 )
            break;
        end
        println( v_i, "..", v_j, " :   Time: ", t.time );
        P = polygon.shortcut( P, v_i, v_j );
    end

    #println( i, "   ", j );
    #println( "Quality: ", quality );

    #shortcut = Polygon2F()
    #push!( shortcut, P[ i ], P[j] );

    #push!( list, P, shortcut );

    bb = BBox2F();
    BBox_bound( bb, P );


    p = BBox_bottom_left( bb );
    for  poly  in list
        Polygon_translate!( poly, p );
    end

    output_polygons_to_file( list, "curves.pdf", true );
end

function  make_straight( P::Polygon{D,T}, quality, ε ) where{D,T}
    if  length( P ) <= 2
        return  P
    end
    if  length( P ) == 3
        q = total_length( P ) / Dist( P[1], P[3] )
        if   ( q <= quality )
            return  P
        end
        Q = Polygon{D,T}() ;
        push!( Q, P[1], P[3] );

        return  Q;
    end

    ℓ = length( P );
    mid::Int = round( Int,  ℓ / 2 );

    PL = cut( P, 1:mid );
    PR = cut( P, mid+1:ℓ )

    QL = make_straight( PL, quality, ε );
    QR = make_straight( PR, quality, ε );
    prefix = length( QL ) ;
    Q = polygon.append_smart!( QL, QR );

    # Now we look for a good enough shortcut between the two parts.
    p_i, q_i, qlt = approx_bichromatic_shortcut( Q, prefix, ε, quality );
    #println( "qlt: ", qlt );
    #println( "shortcut: [",p_i, ":",  q_i, "]" );
    if  qlt > quality
        #println( "Bingo!" );
        Q = polygon.shortcut( Q, p_i, q_i );
    end

    return  Q;
end

function (@main)(ARGS)
    if  length( ARGS ) != 2
        println( "Usage:\n\t"
                 * "shortcut.jl [eps] [file.txt]\n\n" );
        return  -1;
    end

    list = VecPolygon2F();
    ε =  parse(Float64, ARGS[ 1  ] );
    P = read_file( ARGS[ 2 ] );


    println( "N = ", length( P ) );
    push!( list, deepcopy( P ) );

    P_simp_10, _ = frechet_simplify_to_cardin( P, 10 );
    P_simp_30, _ = frechet_simplify_to_cardin( P, 30 );

    Q_800 = make_straight( P, 800.0, ε );
    push!( list, Q_800 );
    Q_400 = make_straight( P, 400.0, ε );
    push!( list, Q_400 );
    Q_200 = make_straight( P, 200.0, ε );
    push!( list, Q_200 );
    Q_100 = make_straight( P, 100.0, ε );
    push!( list, Q_100 );
    Q_50 = make_straight( P, 50.0, ε );
    push!( list, Q_50 );
    Q_25 = make_straight( P, 25.0, ε );
    push!( list, Q_25 );
    Q_10 = make_straight( P, 10.0, ε );
    push!( list, Q_10 );
    Q_2 = make_straight( P, 2.0, ε );
    push!( list, Q_2 );

    output_polygons_to_file( list, "curves_all.pdf", true );

    c,cr,bb = cairo_setup( "curves.pdf", VecPolygon2F( [P] ), true  );

    output_polygons( cr,bb, [P], false, false, 3.0, "Input: "*string(length(P)) );
    Cairo.show_page( cr )

    output_polygons( cr,bb, [P, Q_800], false, false, 3.0, "Straightness: 800" );
    Cairo.show_page( cr )
    output_polygons( cr,bb, [P, Q_400], false, false, 3.0, "Straightness: 400" );
    Cairo.show_page( cr )
    output_polygons( cr,bb, [P, Q_200], false, false, 3.0, "Straightness: 200" );
    Cairo.show_page( cr )
    output_polygons( cr,bb, [P, Q_100], false, false, 3.0, "Straightness: 100" );
    Cairo.show_page( cr )
    output_polygons( cr,bb, [P, Q_50], false, false, 3.0, "Straightness: 50" );
    Cairo.show_page( cr )
    output_polygons( cr,bb, [P, Q_25], false, false, 3.0, "Straightness: 25" );
    Cairo.show_page( cr )
    output_polygons( cr,bb, [P, Q_10], false, false, 3.0, "Straightness: 10" );
    Cairo.show_page( cr )
    output_polygons( cr,bb, [P, Q_2], false, false, 3.0, "Straightness: 2" );
    Cairo.show_page( cr )

    output_polygons( cr,bb, [P, P_simp_10], false, false, 3.0, "Frechet simp 10" );
    Cairo.show_page( cr )
    output_polygons( cr,bb, [P, P_simp_30], false, false, 3.0, "Frechet simp 30" );
    Cairo.show_page( cr )

    Cairo.finish(c);

    return  0;
end


####################################################################3

#export  BBTree_build;

#end

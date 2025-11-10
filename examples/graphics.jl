# Graphics.jl

using Parameters
using Distributions;
using Cairo
using LinearAlgebra
using Printf
using Plots
using PrettyTables
using Dates

using FrechetDist.cg.bbox

VecFloat = Vector{Float64};
VecVecFloat = Vector{VecFloat};


function  draw_disks( cr, disks )
    for  d in disks
        cen = Point2F( d[ 1 ], d[ 2 ] );
        r = d[ 3 ];
        arc( cr, cen[1], cen[2], r, 0.0, 2.0 * pi );
        stroke_preserve(cr);
        Cairo.fill(cr);
    end
end

function  draw_circles( cr, disks )
    for  d in disks
        cen = Point2F( d[ 1 ], d[ 2 ] );
        r = d[ 3 ];
        arc( cr, cen[1], cen[2], r, 0.0, 2.0 * pi );
        
        # Apply the stroke to the path (the circle drawn with arc)
        Cairo.stroke(cr)
    end
end

function  draw_hippodrome( cr, p::Point2F, q::Point2F, r::Float64 )
    r = r * 4.0;

    arc( cr, p[1], p[2], r, 0.0, 2.0 * pi );
    fill_preserve( cr );
    arc( cr, q[1], q[2], r, 0.0, 2.0 * pi );
    fill_preserve( cr );

    v = (p - q) / Dist( p, q );
    u = npoint( -v[2], v[1] );
    px = p - r * u;
    qx = q - r * u;
    qy = q + r * u;
    py = p + r * u;

    move_to(cr, px[ 1 ], px[ 2 ] );
    line_to(cr, qx[ 1 ], qx[ 2 ] );
    line_to(cr, qy[ 1 ], qy[ 2 ] );
    line_to(cr, py[ 1 ], py[ 2 ] );
    line_to(cr, px[ 1 ], px[ 2 ] );
    close_path(cr);
    stroke_preserve(cr);
    Cairo.fill(cr);
    Cairo.stroke( cr );
end

function  draw_polygon_w_offs( cr, P::Polygon2F, offs::VecFloat )
    nv::Int64 = cardin( P );
    for  i in 2:nv
        p = P.pnts[ i - 1 ];
        #p = po.x;
        q = P.pnts[ i ];
        #q = qo.x;

        r = max( offs[ i - 1 ], offs[ i ] );
        #println( " r: ", r );

        get_color_rgbsource_rgb(cr, 1.0, 0.0, 0.8);
        #arc( cr, p[1], p[2], r, 0.0, 2.0 * pi );
        #fill_preserve( cr );

        draw_hippodrome( cr, p, q, r );

        #Cairo.stroke( cr );
        #cr->set_source_rgba(0.0, 0.0, 0.8, 1.0);
        #move_to( cr,  p[1], p[2] )
        #line_to( cr, q[1], q[2] );
        #println( q[1], " ", q[2] );
        Cairo.stroke( cr );
    end
end


function  draw_segments( cr, segs )
    for  s in segs
        p = s.p;
        q = s.q;
        move_to(cr, p[ 1 ], p[ 2 ] );
        line_to(cr, q[ 1 ], q[ 2 ] );
        Cairo.stroke( cr );
    end
end

function  draw_polygon( cr, P, f_close::Bool = false )
    nv::Int64 = cardin( P );
    for  i in 2:nv
        p = P.pnts[ i - 1 ];
        #p = po.x;
        q = P.pnts[ i ];
        #q = qo.x;
        move_to( cr,  p[1], p[2] )
        line_to( cr, q[1], q[2] );
        #println( q[1], " ", q[2] );
    end
    if  (nv > 0 )  &&  ( f_close )
        q = P.pnts[ 1 ];
        line_to( cr, q[1], q[2] );
    end
    Cairo.stroke( cr );
end


function  fill_polygon( cr, P )
    nv::Int64 = cardin( P );
    s = P[1]
    move_to( cr,  s[1], s[2] )
    for  i in 2:nv
        q = P.pnts[ i ];
        line_to( cr, q[1], q[2] );
    end
    if  nv > 0 
        q = P.pnts[ 1 ];
        line_to( cr, q[1], q[2] );
    end
    close_path(cr);
    Cairo.fill(cr)
end


function  draw_points( cr, P::Vector{Point2F}, rad::Float64 = 0.005 )
    nv::Int64 = length( P );
    #rad::Float64 = 0.005
    #println( "nv:", nv );
    for  i in 1:nv
        p = P[ i ];
        move_to( cr, p[1] - rad , p[2] )
        
        arc( cr, p[1], p[2], rad, 0.0, 2 * 3.1415 );
        stroke_preserve(cr);
        fill(cr);
    end
    Cairo.stroke( cr );
end

function  draw_polygon_vertices( cr, P, r::Float64 )
    nv::Int64 = cardin( P );
    #r = total_length( P ) / (100.0*cardin( P ));
    for  i in 1:nv
        p = P.pnts[ i ];
#        set_line_width(cr, 2.00);
#        Cairo.set_source_rgb( cr, 0, 0, 0);
#        println( "RRR = ", r );
        Cairo.arc( cr, p[1], p[2], r, 0.0, 2 * pi);
        Cairo.fill(cr);
#        Cairo.arc( cr, p[1], p[2], r, 0.0, 2 * pi);
#        Cairo.set_source_rgb( cr, 1.0, 1.0, 0);
    end
    Cairo.stroke( cr );
end

function  draw_bbox( cr, bb, scale )
    pa = cg.bbox.bottom_left( bb );
    pc = cg.bbox.top_right( bb );

    pb = npoint( pc[1], pa[2] );
    pd = npoint( pa[1], pc[2] )

    pa = pa * scale;
    pb = pb * scale;
    pc = pc * scale;
    pd = pd * scale;
    ;
    move_to( cr, pa[1], pa[2] )
    line_to( cr, pb[1], pb[2] );
    line_to( cr, pc[1], pc[2] );
    line_to( cr, pd[1], pd[2] );
    line_to( cr, pa[1], pa[2] );

    Cairo.stroke( cr );
end

function  fill_bbox( cr, bb, scale )
    pa = bottom_left( bb );
    pc = top_right( bb );

    pb = npoint( pc[1], pa[2] );
    pd = npoint( pa[1], pc[2] )

    pa = pa * scale;
    pb = pb * scale;
    pc = pc * scale;
    pd = pd * scale;

    move_to( cr, pa[1], pa[2] )
    line_to( cr, pb[1], pb[2] );
    line_to( cr, pc[1], pc[2] );
    line_to( cr, pd[1], pd[2] );
    line_to( cr, pa[1], pa[2] );
    close_path(cr);

    Cairo.fill( cr );
end

function  draw_bboxes( cr, boxes, scale = 1.0 )
    for  bb ∈ boxes
        draw_bbox( cr, bb, scale );
    end
end
function  fill_bboxes( cr, boxes, scale = 1.0 )
    for  bb ∈ boxes
        fill_bbox( cr, bb, scale );
    end
end

function  fill_bboxes_cycle_colors( cr, boxes; α = 1.0, scale = 1.0,
    f_draw_frames = false
)
    count = 0;
    for  bb ∈ boxes
        count += 1
        set_source_rgba( cr, get_color_rgb( count )..., α );

        fill_bbox( cr, bb, scale );
        if  f_draw_frames
            set_source_rgba( cr, get_color_rgb( count )..., 1.0 );
            draw_bbox( cr, bb, scale );
        end
    end
end

function  get_color_rgb( i::Int64 )
    @assert( i > 0 );
    colors =[ 0.0 0.0 0.0;
              1.0 0.0 0.0;
              0.0 1.0 0.0;
              0.0 0.0 1.0;
              0.0 0.5 0.0;
              1.0 0.0 1.0;
              0.0 1.0 1.0;
              1.0 1.0 0.0;
              0.5 0.5 1.0;
              0.5 1.0 0.5;
              0.3 0.5 0.7;
              0.7 0.2 0.4;
              0.5 0.1 0.5 ];
    ind = 1 + ( ( i - 1 ) % size(colors,1));
    #println( ind, " ::: ", size(colors,1) );
    return  colors[ ind, : ];
end


function  compute_bounding_boxes( list::VecPolygon2F )
    bb::BBox2F = BBox2F();

    bound( bb, list );
    expand( bb, 1.05 );
    bbo::BBox2F = deepcopy( bb );
    expand( bbo, 1.05 );

    return  bb, bbo
end

function  get_image_dims( bbo )

    width::Float64 = 1024.0;
    theight::Float64 = 0.0;

    while ( true )
        theight = width * cg.bbox.width( bbo, 2 ) / cg.bbox.width( bbo, 1 );
        if  theight < 2048.0
            break;
        end

        width = width / 2.0;
    end

    iheight::Int64 = convert( Int64, 16 * ceil( theight / 16 ) )
    iwidth::Int64 = convert( Int64, 16 * ceil( width / 16 ) )

    return  iheight,iwidth;
end

function  set_transform( cr, iwidth::Int64, iheight::Int64,
                         bbo::BBox2F )
    xcal = convert( Float64, iwidth) / cg.width( bbo, 1 );

#        ycal = convert( Float64, iheight) / cg.width( bbo, 2 );
#    if  ( ( (xcal / 5.0 ) < ycal ) && ( ycal <= (xcal * 5.0 ) ) )
#        ycal = xcal;
#    end

    Cairo.scale( cr, xcal, xcal );
    bl = cg.bbox.bottom_left( bbo );
    Cairo.translate( cr, -bl[ 1 ], -bl[ 2 ]);
end

function  cairo_setup( filename::String, list::VecPolygon2F,
                       f_pdf::Bool = true )
    bb, bbo = compute_bounding_boxes( list );
    iheight, iwidth = get_image_dims( bbo );

#    set_transform( cr, iwidth, bbo );

    local c
    if (  f_pdf )
        c = Cairo.CairoPDFSurface( filename, iwidth, iheight );
    else
        c = CairoRGBSurface( iwidth, iheight );
    end
    cr = CairoContext(c);

    set_transform( cr, iwidth, iheight, bbo );

    if  ( ! f_pdf )
        set_source_rgb(cr, 1, 1, 1);
        paint(cr);
    end

    return  c,cr,bb;
end

function  cairo_setup( filename::String, top_right::Point2F,
                       f_pdf::Bool = true )
    P = Polygon2F();
    push!( P, Point2F( 0.0, 0.0 ) );
    push!( P, top_right );

    return cairo_setup( filename, [P], f_pdf );
    
end

function  output_polygons( cr,bb::BBox2F,
    list::VecPolygon2F,
    f_draw_vertices::Bool = false,
    f_matching::Bool = false,
    u_width::Float64 = 3.0,
    caption = ""
)
    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    len = length( list );
    count::Int64 = 0;
    for  poly in  list
        count = count + 1;
        #println( count, " ", len );
        set_line_width(cr, u_width );
        #if  len == 2  &&  count == 2
        #    set_source_rgb(cr, 0.0, 0.0, 1.0 );
        #else
        set_source_rgb( cr, get_color_rgb( count )... );

        draw_polygon( cr, poly );
    end

    if  ( f_matching )  &&  ( cardin( list[ 1 ] ) ==  cardin( list[ 2 ] ) )
        P = list[ 1 ];
        Q = list[ 2 ];
        set_line_width(cr, 0.5*u_width);
        set_source_rgb( cr, 1.0, 0.0, 1.0 );
        for  i  in 1:cardin( P )
            p = P[ i ];
            q = Q[ i ];

            move_to( cr, p[1], p[2] )
            line_to( cr, q[1], q[2] );
            Cairo.stroke( cr );
        end
    end

    if  ( f_draw_vertices )
        set_line_width(cr, 0.9*u_width);
        set_source_rgb( cr, 1.0, 0.0, 0.0 );
        for  poly in  list
            draw_polygon_vertices( cr, poly, cg.width( bb) / 200  );
        end
    end

    if  length( caption ) > 0
        # Set font
        font_face::String = "Sans" #
        Cairo.select_font_face( cr, font_face, Cairo.FONT_SLANT_NORMAL,
            Cairo.FONT_WEIGHT_NORMAL)
        w = cg.width( bb );
        font_size = w/40
        Cairo.set_font_size( cr, font_size )

        # Move to the position where text will be drawn
        pc = cg.bbox.bottom_left( bb );
        Cairo.move_to( cr, pc[1] + font_size, pc[2] + font_size );
        
        # Show the text
        Cairo.set_source_rgb( cr, 0,0.3, 0.3)
        Cairo.show_text(cr, caption)
    end
end


function  output_polygons_to_file(  list::VecPolygon2F, filename,
                                    f_pdf::Bool,
                                    f_draw_vertices::Bool = false,
                                    f_matching::Bool = false;
    u_width::Float64 = 3.0,
    caption = ""
)
    c,cr,bb = cairo_setup( filename, list, f_pdf );


    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    len = length( list );
    count::Int64 = 0;
    for  poly in  list
        count = count + 1;
        #println( count, " ", len );
        set_line_width(cr, u_width );
        #if  len == 2  &&  count == 2
        #    set_source_rgb(cr, 0.0, 0.0, 1.0 );
        #else
        set_source_rgb( cr, get_color_rgb( count )... );

        draw_polygon( cr, poly );
    end

    if  ( f_matching )  &&  ( cardin( list[ 1 ] ) ==  cardin( list[ 2 ] ) )
        P = list[ 1 ];
        Q = list[ 2 ];
        set_line_width(cr, 0.5*u_width);
        set_source_rgb( cr, 1.0, 0.0, 1.0 );
        for  i  in 1:cardin( P )
            p = P[ i ];
            q = Q[ i ];

            move_to( cr, p[1], p[2] )
            line_to( cr, q[1], q[2] );
            Cairo.stroke( cr );
        end
    end

    if  ( f_draw_vertices )
        set_line_width(cr, 0.9*u_width);
        set_source_rgb( cr, 1.0, 0.0, 0.0 );
        for  poly in  list
            draw_polygon_vertices( cr, poly, cg.width( bb) / 200  );
        end
    end

    if  length( caption ) > 0
        # Set font
        Cairo.select_font_face(ctx, font_face, Cairo.FONT_SLANT_NORMAL,
            Cairo.FONT_WEIGHT_NORMAL)
        Cairo.set_font_size(ctx, 16 )

        # Move to the position where text will be drawn
        Cairo.move_to( cr, 0, 20)
        
        # Show the text
        Cairo.show_text(cr, caption)
    end

    if  ( ! f_pdf )
        Cairo.write_to_png( c, filename );
    end
    Cairo.finish(c);
end


function  output_polygons_to_file_with_offsets(
    list::VecPolygon2F,
    loffs::VecVecFloat,
    filename,
    f_pdf::Bool,
    f_draw_vertices::Bool = false,
    f_matching::Bool = false
    )

    c,cr,bb = cairo_setup( filename, list, f_pdf );

    u_width::Float64 = 1024.0 * (cg.width( bb) / 4500.0);

    cg.bbox.print( bb );
    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
#    set_line_width(cr, 10.0);
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    len = length( list );
    count::Int64 = 0;
    off::VecFloat = VecFloat();
    for  i in  eachindex(list)
        poly = list[ i ];
        f_off::Bool = false;
        if  ( i <= length( loffs ) )
            f_off = true;
            off = loffs[ i ];
        end
        count = count + 1;
        #println( count, " ", len );
        set_line_width(cr, u_width );
        if  len == 2  &&  count == 2
            set_source_rgb(cr, 0.0, 0.0, 1.0 );
        else
            set_source_rgb( cr, 0.0, 1.0, 0.0 );
        end

        if  ( f_off )
            draw_polygon_w_offs( cr, poly, off );
        else
            draw_polygon( cr, poly );
        end
    end

    if  ( f_matching )  &&  ( cardin( list[ 1 ] ) ==  cardin( list[ 2 ] ) )
        P = list[ 1 ];
        Q = list[ 2 ];
        set_line_width(cr, 1.5*u_width);
        set_source_rgb( cr, 1.0, 0.0, 1.0 );
        for  i  in 1:cardin( P )
            p = P[ i ];
            q = Q[ i ];

            move_to( cr, p[1], p[2] )
            line_to( cr, q[1], q[2] );
            Cairo.stroke( cr );
        end
    end

    if  ( f_draw_vertices )
        set_line_width(cr, 2.0*u_width);
        set_source_rgb( cr, 1.0, 0.0, 0.0 );
        for  poly in  list
            draw_polygon_vertices( cr, poly, cg.width( bb) / 200  );
        end
    end


    if  ( ! f_pdf )
        Cairo.write_to_png( c, filename );
    end
    Cairo.finish(c);
end





#----------------------------------------------------------------
# Output the morphing to a pdf file
function  output_morphing( m::Morphing{N,T}, filename )  where {N,T}
    c,cr,bb = cairo_setup( filename, [ m.P, m.Q ], true );

    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
    set_line_width(cr, 1.0);
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    P = m.P;
    Q = m.Q;
    PS, QS = Morphing_as_polygons( m );

    set_line_width(cr, 0.5);
    set_source_rgb(cr, 1.0, 0.0, 0.0 );

    nv::Int64 = cardin( PS );
    for  i in 1:nv
        #loc::Point2I  = sol.pnts[ i ];

        po::Point{N,T} = PS[ i ]
        qo::Point{N,T} = QS[ i ]

#        println( loc.x[ 1 ], ", ", loc.x[ 2] );

        move_to( cr,  po[1], po[2] )
        line_to( cr, qo[1], qo[2] );
        Cairo.stroke( cr );
    end

    set_line_width(cr, 1.0);
    set_source_rgb(cr, 0.0, 1.0, 0.0 );
    draw_polygon( cr, P );
    set_source_rgb(cr, 0.0, 0.0, 1.0 );
    draw_polygon( cr, Q );

    set_source_rgb( cr, 1,0,0 );

    Cairo.finish(c);
end


function  draw_frames( cr, sp::Segment2F, sq::Segment2F,
                       frames::Int64, P::Polygon2F, Q::Polygon2F, bb::BBox2F )

    delta::Float64 = 1.0 / (frames -1 );
    t::Float64 = 0.0;

    for i in 1:frames
        set_line_width(cr, 3.5);
        set_source_rgb(cr, 0.0, 0.8, 0.0 );
        draw_polygon( cr, P );
        set_source_rgb(cr, 0.0, 0.0, 1.0 );
        draw_polygon( cr, Q );


        set_line_width(cr, 10.5);
        set_source_rgb( cr, 1,0,0 );
        p::Point2F = segment.at( sp, t );
        q::Point2F = segment.at( sq, t );

        move_to( cr,  p[1], p[2]  )
        line_to( cr, q[1], q[2] );
        Cairo.stroke( cr );

        Cairo.show_page( cr );

        t = t + delta;
    end
end


function  compute_frames( pout, qout, total_frames )
#    println( "TOTAL FRAMES :", total_frames );
    lpout::Vector{Float64} = Polygon_prefix_lengths( pout )
    lqout::Vector{Float64} = Polygon_prefix_lengths( qout )

    lens::Vector{Float64} = lpout + lqout;

    steps::Int64 = length( lens ) - 1;
    total_length = last( lens );


    # compute how many frames for each leg...
    frames::Vector{Int64} = zeros( Int64, steps );
    # compute how many frames for each leg...
#    frames::Vector{Int64} = zeros( Int64, steps );
    acc = 0;
    for  i  in  1:steps
        leg_len = lens[ i + 1 ] - lens[ i ];
        if  ( leg_len <= 0 )
            continue;
        end
        acc += leg_len;
        num_pnts = round( Int64, total_frames * leg_len / total_length );
#        println( "num_pnts: ", num_pnts );
        if  num_pnts > 0
            acc = 0;
            frames[ i ] = max( num_pnts, 2 );
            continue;
        end

        # Very short edge...
        num_pnts = round( Int64, total_frames * acc / total_length );
        if  ( num_pnts == 0 )
            continue; # skip it....
        end
        acc = 0;
        frames[ i ] = num_pnts;
    end

    #println( "Total frames #: ", sum( frames ) );
    return  frames, steps
end


function  output_frechet_movie( m::Morphing{N,T},
                                filename::String,
    total_frames::Int64 = 800,
    f_show_vertices::Bool = false ) where {N,T}

    if  isfile( filename )
        println( "\n\n",  filename, " already exists...\n\n\n" );
        return;
    end

    c,cr,bb = cairo_setup( filename, [ m.P, m.Q ], true );

    #   Cairo.save( cr );

    set_line_width(cr, 0.5);
    set_source_rgb(cr, 1.0, 0.0, 0.0 );


    pout,qout = Morphing_as_polygons( m );

    np::Int64 = cardin( pout );
    nq::Int64 = cardin( qout );
    if   np != nq
        println( "Error np!= nq" );
        exit( -1 );
    end

    lpout::Vector{Float64} = Polygon_prefix_lengths( pout )
    lqout::Vector{Float64} = Polygon_prefix_lengths( qout )

    lens::Vector{Float64} = lpout + lqout;

    steps = length( lens ) - 1;
    total_length = last( lens );
    #println( lens );


    # compute how many frames for each leg...
    frames, steps = compute_frames( pout, qout, total_frames )

    total_output_frames = sum( frames );
    println( "Total number of frames: ", total_output_frames );
    skip::Int64 = max( 1, floor( total_output_frames / total_frames ) );
    f_do_skip::Bool = true ;

    for  i in 1:steps
        count = 0;
        if  frames[ i ] == 0
            continue;
        end
        if  ( ! f_do_skip )
            f_output_frame = true;
        else
            f_output_frame = ( count % skip == 0 )  ||  ( i == steps );
        end
        count = count + 1;
        if  ( f_output_frame )
            sp = Segment( pout[ i ], pout[ i + 1 ] );
            sq = Segment( qout[ i ], qout[ i + 1 ] );

            draw_frames( cr, sp, sq, frames[ i ], m.P, m.Q, bb );
        end
    end

    Cairo.finish(c);
end




mutable struct ContextMovie
    bb::BBox2F;
    bbo::BBox2F;
    iheight::Int64
    iwidth::Int64;
    frame_count::Int64;
    dir::String;
    f_show_vertices::Bool
end


mutable struct RecFrame
    frame::Int64
    p::Point2F
    q::Point2F
end
function   draw_image_record( cm::ContextMovie, P, Q, p, q,
                              vec::Vector{RecFrame} )
    cm.frame_count += 1;
    push!( vec, RecFrame( cm.frame_count, deepcopy( p ), deepcopy( q  ) ) );
end


"""
    draw_image_frame
"""
function   draw_image_frame( cm::ContextMovie, P, Q, rf::RecFrame )

    filename = cm.dir * "/" * @sprintf( "%06d.png", rf.frame )

    c = CairoRGBSurface(cm.iwidth, cm.iheight );
    cr = CairoContext(c);


    set_source_rgb(cr, 1, 1, 1);
    paint(cr);

    set_transform( cr, cm.iwidth, cm.iheight, cm.bbo );

    set_line_width(cr, 3.5);
    set_source_rgb(cr, 0.0, 0.8, 0.0 );
    draw_polygon( cr, P );
    if   ( cm.f_show_vertices )
        set_source_rgb(cr, 0.0, 1.0, 0.0 );
        draw_polygon_vertices( cr, P, cg.width( cm.bb ) / 200 );
    end
    set_source_rgb(cr, 0.0, 0.0, 0.8 );
    draw_polygon( cr, Q );
    if   ( cm.f_show_vertices )
        set_source_rgb(cr, 0.0, 0.0, 1.0 );
        draw_polygon_vertices( cr, Q, cg.width( cm.bb ) / 200 );
    end

    set_line_width(cr, 10.5);
    set_source_rgb( cr, 1,0,0 );

    move_to( cr,  rf.p[1], rf.p[2]  )
    line_to( cr, rf.q[1], rf.q[2] );
    Cairo.stroke( cr );

    #print( "filename :", filename, "\r" );
    Cairo.write_to_png( c, filename );
    #Cairo.show_page( cr );
end



function  draw_frames_images( cm::ContextMovie, sp::Segment2F, sq::Segment2F,
    frames::Int64, P::Polygon2F, Q::Polygon2F, vec::Vector{RecFrame}  )

    delta::Float64 = 1.0 / (frames -1 );
    t::Float64 = 0.0;

    for i in 1:frames
        p::Point2F = segment.at( sp, t );
        q::Point2F = segment.at( sq, t );

        draw_image_record( cm, P, Q, p, q, vec )
#        cm::ContextMovie, P, Q, p, q,
#                               )

        #draw_image_frame( cm, P, Q, p, q );

        t = t + delta;
    end
end

function  frames_generate( cm::ContextMovie, P, Q, vec_rf )
    for rf in vec_rf
        draw_image_frame( cm::ContextMovie, P, Q, rf )
    end
    return  1
end



function  rmx( tmp_filename )
    if isfile( tmp_filename )
        rm( tmp_filename );
    end
end


function  output_frechet_movie_mp4( m::Morphing{N,T},
                                filename::String,
    total_frames::Int64 = 800,
    f_show_vertices::Bool = false
) where {N,T}
    f_debug::Bool = false;
    cm = ContextMovie(BBox2F(), BBox2F(), 0, 0, 0, "/tmp/r_draw/",
        f_show_vertices );
    cm.bb, cm.bbo = compute_bounding_boxes( [ m.P, m.Q ] );
    cm.iheight, cm.iwidth = get_image_dims( cm.bbo );

    pout,qout = Morphing_as_polygons( m );

    np::Int64 = cardin( pout );
    nq::Int64 = cardin( qout );
    if   np != nq
        println( "Error np!= nq" );
        exit( -1 );
    end

    #set_transform( cr, iwidth, bbo );

    #set_transform( cr, iwidth, bbo );

    #c = Cairo.CairoPDFSurface( filename, iwidth, iheight );
    #cr = CairoContext(c);
    # Create temporary directory for images...
    if  isdir( cm.dir )
        rm( cm.dir, recursive=true)
    end
    mkdir( cm.dir );


    frames, steps = compute_frames( pout, qout, total_frames )
    total_output_frames = sum( frames );
    skip::Int64 = max( 1, floor( total_output_frames / total_frames ) );
    f_do_skip::Bool = true ;

    vec_rf = Vector{RecFrame}();
    frame_count::Int64 = 0;

    # We first calculate the frames we need... into vec_rf...
    for  i in 1:steps
        count = 0;
        if  frames[ i ] == 0
            continue;
        end
        if  ( ! f_do_skip )
            f_output_frame = true;
        else
            f_output_frame = ( count % skip == 0 )  ||  ( i == steps );
        end
        count = count + 1;

        if  ( f_output_frame )
            sp = Segment( pout[ i ], pout[ i + 1 ] );
            sq = Segment( qout[ i ], qout[ i + 1 ] );

            draw_frames_images( cm, sp, sq, frames[ i ], m.P, m.Q, vec_rf );
        end
    end


    f_debug  &&  println( "Splitting to threads... Sit back ;)" );
    f_debug  &&  println( Threads.nthreads() );

    chunks = Iterators.partition(vec_rf, length(vec_rf) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn frames_generate( cm, m.P, m.Q, chunk );
    end
    chunk_sums = fetch.(tasks)

    ####
    # We now call ffmpeg to create the movie...
    # ffmpeg -r 10 -i temp/ do.mp4

#    cmd = ( " -r 10 -i " *   "
#          * filename );
    tmp_filename = "tmp.mp4";
    rmx( filename );
    rmx( tmp_filename );
#    println( "ffmpeg $options" );
    f_debug && println( "Encoding movie with ffmpeg..." );
    options = [ "-r", "10", "-i",
               cm.dir * "/" * "%06d.png",
               "-c:v", "libx264", tmp_filename ];
    output = read(pipeline( `ffmpeg $options`, stderr="/tmp/errs.txt" ),
        String);

    f_debug &&  println( "Rencoding with handbrake..." );

    # HandBrakeCLI -Z  -i movie.mp4  -o movie_2.mp4
    options_2 = [ "-Z", "Android 1080p30", "-i", tmp_filename,
                 "-o", filename ];
    output_2 = read(pipeline( `HandBrakeCLI $options_2`,
                            stderr="/tmp/errs_2.txt" ), String);
    rmx( tmp_filename );
    if  isfile( filename )
        println( "Created: ", filename );
    end
end


mutable struct ORecFrame
    frame::Int64
    t::Float64
    R::Polygon2F
end


function  frames_generate_o( cm::ContextMovie, P, Q, vec_rf )
    for rf in vec_rf
        draw_image_frame_o( cm::ContextMovie, P, Q, rf )
    end
    return  1
end


function   encode_to_mp4( dir, filename )
    tmp_filename = "tmp.mp4";
    rmx( filename );
    rmx( tmp_filename );
    #println( "Encoding movie with ffmpeg..." );
    options = [ "-r", "10", "-i",
               dir * "/" * "%06d.png",
               "-c:v", "libx264", tmp_filename ];
    output = read(pipeline( `ffmpeg $options`, stderr="/tmp/errs.txt" ),
        String);


    #println( "\n\n\n\n\n\n" );
    #println( "Rencoding with handbrake..." );
    #println( filename );
    #println( "\n\n\n\n\n\n" );

    # HandBrakeCLI -Z  -i movie.mp4  -o movie_2.mp4
    options_2 = [ "-Z", "Android 1080p30", "-i", tmp_filename,
                 "-o", filename ];
    output_2 = read(pipeline( `HandBrakeCLI $options_2`,
                            stderr="/tmp/errs_2.txt" ), String);
    rmx( tmp_filename );
    if  isfile( filename )
        println( "Created: ", filename );
    end
end

mutable struct ORecFrame
    frame::Int64
    t::Float64
    R::Polygon2F
end

"""
    draw_image_frame_o
"""
function   draw_image_frame_o( cm::ContextMovie, P, Q, rf::ORecFrame )

    filename = cm.dir * "/" * @sprintf( "%06d.png", rf.frame )

    c = CairoRGBSurface(cm.iwidth, cm.iheight );
    cr = CairoContext(c);


    set_source_rgb(cr, 1, 1, 1);
    paint(cr);

    set_transform( cr, cm.iwidth, cm.iheight, cm.bbo );

    set_line_width(cr, 1);
    set_source_rgb(cr, 0.3, 0.3, 0.3 );
    draw_polygon( cr, P );
    draw_polygon( cr, Q );


    set_line_width(cr, 4);
    set_source_rgb(cr, 0.0, 0.8, 0.0 );
    draw_polygon( cr, rf.R );
    if   ( cm.f_show_vertices )
        set_source_rgb(cr, 0.0, 1.0, 0.0 );
        draw_polygon_vertices( cr, rf.R, cg.width( cm.bb ) / 200 );
    end
    Cairo.stroke( cr );

    Cairo.write_to_png( c, filename );
end


function  frames_generate_o( cm::ContextMovie, P, Q,
                             vec_rf::Vector{ORecFrame} )
    for rf in vec_rf
        draw_image_frame_o( cm::ContextMovie, P, Q, rf )
    end
    return  1
end




function  output_ortho_frechet_movie_mp4( m::Morphing{N,T},
                                          filename::String,
                                          total_frames::Int64 = 200,
                                          ) where {N,T}
    cm = ContextMovie(BBox2F(), BBox2F(), 0, 0, 0, "/tmp/r_draw/",
                      false  );
    cm.bb, cm.bbo = compute_bounding_boxes( [ m.P, m.Q ] );
    cm.iheight, cm.iwidth = get_image_dims( cm.bbo );

    pout,qout = Morphing_as_polygons( m );

    np::Int64 = cardin( pout );
    nq::Int64 = cardin( qout );
    @assert( np == nq );

    if  isdir( cm.dir )
        rm( cm.dir, recursive=true)
    end
    mkdir( cm.dir );

    vec_rf = Vector{ORecFrame}();

    delta = 1.0 / ( total_frames - 1 );
    t::Float64 = 1.0;

    # We first calculate the frames we need... into vec_rf...
    for  i in 1:total_frames
        R = Polygon_convex_comb( pout, qout, t );

        ofr = ORecFrame( i, t, R );
        push!( vec_rf, ofr );
        t = t - delta;
    end

    #println( "Splitting to threads... Sit back ;)" );

    #println( Threads.nthreads() );

    chunks = Iterators.partition(vec_rf, length(vec_rf) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn frames_generate_o( cm, m.P, m.Q, chunk );
    end
    chunk_sums = fetch.(tasks)

    encode_to_mp4( cm.dir, filename );
end


function  get_diagram_locs( PE::Vector{EventPoint{N,T}}, P::Polygon{N,T}
                            ) where {N,T}

    prefixes::Vector{Float64} = Polygon_prefix_lengths( P )

    pout = Vector{T}()

    len = length( PE );
    push!( pout, 0 );

    i = 2;
    while  ( i <= len )
        ep = PE[ i ];
        if  ep.type == PT_VERTEX
            push!( pout, prefixes[ ep.i ] );
            i = i + 1;
            continue;
        end

        loc = ep.i;
        edge_length = prefixes[ loc + 1 ] - prefixes[ loc ];
        push!( pout, prefixes[ loc ] + ep.t * edge_length );

        i = i + 1
    end

    return  pout
end

#P::Polygon{N,T}, Q::Polygon{N,T}, Pe, Qe,
function  output_frechet_diagram( m::Morphing{N,T}, filename )  where {N,T}

    P_coords::Vector{T} = get_diagram_locs( m.pes, m.P );
    Q_coords::Vector{T} = get_diagram_locs( m.qes, m.Q );

    len = length( P_coords );
    poly = Polygon2F();
    for  i in 1:len
        push_smart!( poly, npoint( P_coords[ i ], Q_coords[ i ] ) );
    end

    psum::Vector{Float64} = Polygon_prefix_lengths( m.P )
    qsum::Vector{Float64} = Polygon_prefix_lengths( m.Q )

    c,cr,bb = cairo_setup( filename, [ poly ], true );

    set_source_rgb(cr,0.9,0.0,0.0);
    set_line_width(cr, 1.0);
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    plen = cardin( m.P );
    qlen = cardin( m.Q );
    ymax = last( qsum );
    xmax = last( psum );

    for  i in 1:plen
        move_to( cr, psum[i], 0 )
        line_to( cr, psum[i], ymax );
        Cairo.stroke( cr );
    end

    for  i in 1:qlen
        move_to( cr,    0, qsum[i] );
        line_to( cr, xmax, qsum[i] );
        Cairo.stroke( cr );
    end
    set_line_width(cr, 1.00);
    set_source_rgb(cr, 0.0, 1.0, 0.0 );
    draw_polygon( cr, poly );
    Cairo.finish(c);
end


function flip_y_axis(ctx::CairoContext, height::Real)
    # Scale the y-axis by -1.
    # This flips the y-axis, but also mirrors the content.
    Cairo.scale(ctx, 1.0, -1.0)
    
    # Translate the context back.
    # The negative height is needed to move the drawing back into view.
    Cairo.translate(ctx, -0.0, -height + 0.01)
end

function add_to_pdf_text( cr, n, text_to_write::String )
    Cairo.save( cr )
    flip_y_axis( cr, n+1 )

    scale = n * 0.02 / 8.0
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


# ????
# End of the frechet distance computation part
#################################################################
#################################################################
#################################################################

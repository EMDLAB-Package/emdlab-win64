% stator tooth & coil

function emdlab_g2d_lib_tc105(g, wst, wss, dsy, dss, frr, xShift, yShift, mirrorFlag, name1, name2)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    wst = 10;
    wss = 12;
    dsy = 7;
    dss = 20;
    frr = 0.7;
    xShift = 0;
    yShift = 0;
    mirrorFlag = 0;
    name1 = 'stator';
    name2 = 'sca';
end

% check fillet ratio
if (frr <= 0)
    emdlab_g2d_lib_tc101(g, wst, wss, dsy, dss, xShift, yShift, mirrorFlag, name1, name2);
    return;
end

if (frr >= 1)

    x1 = 0;
    y1 = 0;

    x2 = wss/2+wst/2;
    y2 = 0;

    x3 = x2;
    y3 = dsy;

    x4 = wst/2;
    y4 = y3 + wss/2;

    x5 = x4;
    y5 = y3 + dss;

    x6 = 0;
    y6 = y5;

    x7 = x2;
    y7 = y5;

    x8 = x2;
    y8 = y4;

    if mirrorFlag
        y1 = -y1;
        y2 = -y2;
        y3 = -y3;
        y4 = -y4;
        y5 = -y5;
        y6 = -y6;
        y7 = -y7;
        y8 = -y8;
    end

    p1 = g.addPoint(x1 + xShift,y1 + yShift);
    p2 = g.addPoint(x2 + xShift,y2 + yShift);
    p3 = g.addPoint(x3 + xShift,y3 + yShift);
    p4 = g.addPoint(x4 + xShift,y4 + yShift);
    p5 = g.addPoint(x5 + xShift,y5 + yShift);
    p6 = g.addPoint(x6 + xShift,y6 + yShift);
    p7 = g.addPoint(x7 + xShift,y7 + yShift);
    p8 = g.addPoint(x8 + xShift,y8 + yShift);

    e1 = g.addSegment(p1,p2);
    e2 = g.addSegment(p2,p3);
    e3 = g.addArc(p8,p3,p4,mirrorFlag);
    e4 = g.addSegment(p4,p5);
    e5 = g.addSegment(p5,p6);
    e6 = g.addSegment(p6,p1);
    e7 = g.addSegment(p3,p7);
    e8 = g.addSegment(p7,p5);

    if mirrorFlag
        l1 = g.addLoop(-e6,-e5,-e4,-e3,-e2,-e1);
        l2 = g.addLoop(e3,e4,-e8,-e7);
    else
        l1 = g.addLoop(e1,e2,e3,e4,e5,e6);
        l2 = g.addLoop(e7,e8,-e4,-e3);
    end

else

    x1 = 0;
    y1 = 0;

    x2 = wss/2+wst/2;
    y2 = 0;

    x3 = x2;
    y3 = dsy;

    % fillet radius
    r_fillet = frr * wss/2;

    x4 = wst/2 + r_fillet;
    y4 = y3;

    x5 = wst/2;
    y5 = y3 + r_fillet;

    x6 = x5;
    y6 = dsy + dss;

    x7 = 0;
    y7 = y6;

    x8 = x2;
    y8 = y6;

    x9 = x4;
    y9 = y5;

    if mirrorFlag
        y1 = -y1;
        y2 = -y2;
        y3 = -y3;
        y4 = -y4;
        y5 = -y5;
        y6 = -y6;
        y7 = -y7;
        y8 = -y8;
        y9 = -y9;
    end

    p1 = g.addPoint(x1 + xShift,y1 + yShift);
    p2 = g.addPoint(x2 + xShift,y2 + yShift);
    p3 = g.addPoint(x3 + xShift,y3 + yShift);
    p4 = g.addPoint(x4 + xShift,y4 + yShift);
    p5 = g.addPoint(x5 + xShift,y5 + yShift);
    p6 = g.addPoint(x6 + xShift,y6 + yShift);
    p7 = g.addPoint(x7 + xShift,y7 + yShift);
    p8 = g.addPoint(x8 + xShift,y8 + yShift);
    p9 = g.addPoint(x9 + xShift,y9 + yShift);

    e1 = g.addSegment(p1,p2);
    e2 = g.addSegment(p2,p3);
    e3 = g.addSegment(p3,p4);
    e4 = g.addArc(p9,p4,p5,mirrorFlag);
    e5 = g.addSegment(p5,p6);
    e6 = g.addSegment(p6,p7);
    e7 = g.addSegment(p7,p1);
    e8 = g.addSegment(p3,p8);
    e9 = g.addSegment(p8,p6);

    if mirrorFlag
        l1 = g.addLoop(-e7,-e6,-e5,-e4,-e3,-e2,-e1);
        l2 = g.addLoop(e3,e4,e5,-e9,-e8);
    else
        l1 = g.addLoop(e1,e2,e3,e4,e5,e6,e7);
        l2 = g.addLoop(e8,e9,-e5,-e4,-e3);
    end

end

g.addFace(name1, l1);
g.addFace(name2, l2);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,255,137,39);

% visualizations for debug
close all;
g.setMeshMaxLength(0.2);
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
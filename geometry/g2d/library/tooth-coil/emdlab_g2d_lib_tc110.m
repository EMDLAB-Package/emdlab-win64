% stator tooth & coil

function emdlab_g2d_lib_tc110(g, wst, wss, dss, bs0, hs0, tta, xShift, yShift, mirrorFlag, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    wst = 10;
    wss = 18;
    dss = 40;
    bs0 = 5;
    hs0 = 2;
    tta = 15;
    xShift = 0;
    yShift = 0;
    mirrorFlag = false;
    name1 = 'stator';
    name2 = 'sca';
    name3 = 'sap';
end

tta = max(tta,0);
tta = deg2rad(tta);
bs0 = min(bs0,0.95*wss);
bs0 = max(bs0,0.05*wss);
ltmp = wss/2-bs0/2;
hs1 = ltmp*tan(tta) + hs0;

x1 = 0;
y1 = 0;

x2 = wst/2;
y2 = 0;

x3 = x2;
y3 = dss/2 - hs1;

x4 = wss/2 + wst/2 - bs0/2;
y4 = dss/2 - hs0;

x5 = x4;
y5 = dss/2;

x6 = 0;
y6 = y5;

x7 = wss/2 + wst/2;
y7 = 0;

x8 = x7;
y8 = y3;

x9 = x7;
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
e4 = g.addSegment(p4,p5);
e5 = g.addSegment(p5,p6);
e6 = g.addSegment(p6,p1);
e7 = g.addSegment(p2,p7);
e8 = g.addSegment(p7,p8);
e9 = g.addSegment(p8,p9);
e10 = g.addSegment(p9,p5);

if mirrorFlag

    l1 = g.addLoop(-e6,-e5,-e4,-e3,-e2,-e1);
    if tta == 0
        e11 = g.addSegment(p8,p4);
        l2 = g.addLoop(e2,e3,-e11,-e8,-e7);
        l3 = g.addLoop(e11,e4,-e10,-e9);
    else
        e11 = g.addSegment(p8,p3);
        l2 = g.addLoop(e2,-e11,-e8,-e7);
        l3 = g.addLoop(e11,e3,e4,-e10,-e9);
    end

else

    l1 = g.addLoop(e1,e2,e3,e4,e5,e6);
    if tta == 0
        e11 = g.addSegment(p8,p4);
        l2 = g.addLoop(e7,e8,e11,-e3,-e2);
        l3 = g.addLoop(e9,e10,-e4,-e11);
    else
        e11 = g.addSegment(p8,p3);
        l2 = g.addLoop(e7,e8,e11,-e2);
        l3 = g.addLoop(e9,e10,-e4,-e3,-e11);
    end

end

g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3, l3);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,255,137,39);
g.setFaceColor(name3,0,255,255);

% visualizations for debug
close all;
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
% stator tooth & coil

function emdlab_g2d_lib_tc51(g, wst, wss, dsy, dss, bs0, hs0, tta, xShift, yShift, mirrorFlag, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    wst = 10;
    wss = 18;
    dss = 20;
    dsy = 10;
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

x2 = wss/2 + wst/2;
y2 = 0;

x3 = x2;
y3 = dsy;

x4 = wst/2;
y4 = dsy;

x5 = x4;
y5 = dsy + dss - hs1;

x6 = wst/2 + wss/2 - bs0/2;
y6 = dsy + dss - hs0;

x7 = x6;
y7 = y6 + hs0;

x8 = 0;
y8 = y7;

x9 = x2;
y9 = y5;

x10 = x2;
y10 = y7;

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
    y10 = -y10;
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
p10 = g.addPoint(x10 + xShift,y10 + yShift);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addSegment(p3,p4);
e4 = g.addSegment(p4,p5);
e5 = g.addSegment(p5,p6);
e6 = g.addSegment(p6,p7);
e7 = g.addSegment(p7,p8);
e8 = g.addSegment(p8,p1);
e9 = g.addSegment(p3,p9);
e10 = g.addSegment(p9,p10);
e11 = g.addSegment(p10,p7);

if mirrorFlag

    l1 = g.addLoop(-e8,-e7,-e6,-e5,-e4,-e3,-e2,-e1);
    if tta == 0
        e12 = g.addSegment(p9,p6);
        l2 = g.addLoop(e4,e5,-e12,-e9,e3);
        l3 = g.addLoop(e12,e6,-e11,-e10);
    else
        e12 = g.addSegment(p9,p5);
        l2 = g.addLoop(e3,e4,-e12,-e9);
        l3 = g.addLoop(e12,e5,e6,-e11,-e10);
    end

else

    l1 = g.addLoop(e1,e2,e3,e4,e5,e6,e7,e8);
    if tta == 0
        e12 = g.addSegment(p9,p6);
        l2 = g.addLoop(-e3,e9,e12,-e5,-e4);
        l3 = g.addLoop(e10,e11,-e6,-e12);
    else
        e12 = g.addSegment(p9,p5);
        l2 = g.addLoop(e9,e12,-e4,-e3);
        l3 = g.addLoop(e10,e11,-e6,-e5,-e12);
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
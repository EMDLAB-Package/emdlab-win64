% stator tooth & coil

function emdlab_g2d_lib_tc109(g, wst, wss, dss, xShift, yShift, mirrorFlag, name1, name2)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    wst = 10;
    wss = 18;
    dss = 20;
    xShift = 0;
    yShift = 0;
    mirrorFlag = false;
    name1 = 'stator';
    name2 = 'sca';
end

x1 = 0;
y1 = 0;

x2 = wst/2;
y2 = 0;

x3 = x2;
y3 = dss/2;

x4 = 0;
y4 = y3;

x5 = wst/2+wss/2;
y5 = 0;

x6 = x5;
y6 = y4;

if mirrorFlag
    y1 = -y1;
    y2 = -y2;
    y3 = -y3;
    y4 = -y4;
    y5 = -y5;
    y6 = -y6;
end

p1 = g.addPoint(x1 + xShift,y1 + yShift);
p2 = g.addPoint(x2 + xShift,y2 + yShift);
p3 = g.addPoint(x3 + xShift,y3 + yShift);
p4 = g.addPoint(x4 + xShift,y4 + yShift);
p5 = g.addPoint(x5 + xShift,y5 + yShift);
p6 = g.addPoint(x6 + xShift,y6 + yShift);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addSegment(p3,p4);
e4 = g.addSegment(p4,p1);
e5 = g.addSegment(p2,p5);
e6 = g.addSegment(p5,p6);
e7 = g.addSegment(p6,p3);

if mirrorFlag
    l1 = g.addLoop(-e4,-e3,-e2,-e1);
    l2 = g.addLoop(e2,-e7,-e6,-e5);
else
    l1 = g.addLoop(e1,e2,e3,e4);
    l2 = g.addLoop(e5,e6,e7,-e2);
end

g.addFace(name1, l1);
g.addFace(name2, l2);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,255,137,39);

% visualizations for debug
close all;
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
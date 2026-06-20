% rotor & magnet
% inner rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_rm_spm102(g, wpole, dry, dm, embrace, xShift, yShift, mirrorFlag, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    wpole = 20;
    dry = 5;
    dm = 4;
    embrace = 0.8;    
    xShift = 0;
    yShift = 0;
    mirrorFlag = false;
    name1 = 'rotor';
    name2 = 'magnet';
    name3 = 'rap';
end

x1 = 0;
y1 = 0;

x2 = wpole/2;
y2 = 0;

x3 = x2;
y3 = dry;

x4 = embrace*wpole/2;
y4 = y3;

x5 = 0;
y5 = y4;

x6 = x2;
y6 = dm + dry;

x7 = x4;
y7 = y6;

x8 = 0;
y8 = y6;

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
e3 = g.addSegment(p3,p4);
e4 = g.addSegment(p4,p5);
e5 = g.addSegment(p5,p1);
e6 = g.addSegment(p3,p6);
e7 = g.addSegment(p6,p7);
e8 = g.addSegment(p7,p8);
e9 = g.addSegment(p8,p5);
e10 = g.addSegment(p4,p7);

if mirrorFlag

    l1 = g.addLoop(-e5,-e4,-e3,-e2,-e1);
    l2 = g.addLoop(e4,-e9,-e8,-e10);
    l3 = g.addLoop(e3,e10,-e7,-e6);

else

    l1 = g.addLoop(e1,e2,e3,e4,e5);
    l2 = g.addLoop(e10,e8,e9,-e4);
    l3 = g.addLoop(e6,e7,-e10,-e3);

end

g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3, l3);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,160,78,146);
g.setFaceColor(name3,0,255,255);

% visualizations for debug
close all;
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
% rotor & magnet
% outer rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_rm_ipm2(g, ID, OD, p, dm, wtrib, wrrib, bm0, t0, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    ID = 50;
    OD = 90;
    p = 8;
    dm = 4;
    wtrib = 1;
    wrrib = 0;
    bm0 = 3;
    t0 = 1;
    name1 = 'rotor';
    name2 = 'magnet';
    name3 = 'rap';
end

% pole pitch angle
alpha_p = 2*pi/p;

x1 = (ID/2);
y1 = 0;

x4 = (OD/2);
y4 = 0;

x5 = (OD/2) * cos(alpha_p/2);
y5 = (OD/2) * sin(alpha_p/2);

x6 = (ID/2) * cos(alpha_p/2);
y6 = (ID/2) * sin(alpha_p/2);

r = OD/2 - wtrib;
[x9,y9] = g.getIntersectionRayCircle(0,-0.5*bm0/cos(alpha_p/2),cos(alpha_p/2),sin(alpha_p/2),0,0,r);
ux = -cos(alpha_p/2);
uy = -sin(alpha_p/2);

err_fcn = @(t) abs(sqrt(r^2 - (y9 + t*uy)^2) - (x9 + t*ux) - dm - t0);
t = fminbnd(err_fcn,0,ID/2);

x7 = (x9 + t*ux);
y7 = (y9 + t*uy);

x8 = x7 + dm;
y8 = y7;

[x10,y10] = g.getIntersectionRayCircle(x8,y8,1,0,0,0,r);

x2 = x7;
y2 = wrrib/2;

x3 = x8;
y3 = wrrib/2;

p1 = g.addPoint(x1,y1);
p2 = g.addPoint(x2,y2);
p3 = g.addPoint(x3,y3);
p4 = g.addPoint(x4,y4);
p5 = g.addPoint(x5,y5);
p6 = g.addPoint(x6,y6);
p7 = g.addPoint(x7,y7);
p8 = g.addPoint(x8,y8);
p9 = g.addPoint(x9,y9);
p10 = g.addPoint(x10,y10);
o = g.addPoint(0,0);

if wrrib == 0

    e1 = g.addSegment(p1,p2);
    e2 = g.addSegment(p2,p3);
    e3 = g.addSegment(p3,p4);
    e4 = g.addArc(o,p4,p5,1);
    e5 = g.addSegment(p5,p6);
    e6 = g.addArc(o,p6,p1,0);
    e7 = g.addSegment(p3,p8);
    e8 = g.addSegment(p8,p7);
    e9 = g.addSegment(p7,p2);
    e10 = g.addArc(o,p10,p9,1);
    e11 = g.addSegment(p9,p7);
    e12 = g.addSegment(p8,p10);

    l1 = g.addLoop(e1,-e9,-e11,-e10,-e12,-e7,e3,e4,e5,e6);
    l2 = g.addLoop(e2,e7,e8,e9);
    l3 = g.addLoop(-e8,e12,e10,e11);

    g.addFace(name1, l1);
    g.addFace(name2, l2);
    g.addFace(name3, l3);

else

    e1 = g.addSegment(p1,p4);
    e2 = g.addArc(o,p4,p5,1);
    e3 = g.addSegment(p5,p6);
    e4 = g.addArc(o,p6,p1,0);
    e5 = g.addSegment(p2,p3);
    e6 = g.addSegment(p3,p8);
    e7 = g.addSegment(p8,p7);
    e8 = g.addSegment(p7,p2);
    e9 = g.addArc(o,p10,p9,1);
    e10 = g.addSegment(p9,p7);
    e11 = g.addSegment(p8,p10);

    l1 = g.addLoop(e1,e2,e3,e4);
    l2 = g.addLoop(e5,e6,e11,e9,e10,e8);
    l3 = g.addLoop(e5,e6,e7,e8);
    l4 = g.addLoop(e9,e10,-e7,e11);

    g.addFace(name1, l1,l2);
    g.addFace(name2, l3);
    g.addFace(name3, l4);

end

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,160,78,146);
g.setFaceColor(name3,0,255,255);

% visualizations for debug
close all;
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
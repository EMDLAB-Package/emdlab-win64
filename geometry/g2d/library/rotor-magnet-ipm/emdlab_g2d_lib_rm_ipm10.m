% rotor & magnet
% outer rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_rm_ipm10(g, ID, OD, p, dm, alpha_v, wtrib, wrrib, bm0, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    ID = 50;
    OD = 90;
    p = 8;
    dm = 4;
    alpha_v = 150;
    wtrib = 1;
    wrrib = 1;
    bm0 = 3;
    name1 = 'rotor';
    name2 = 'magnet';
    name3 = 'rap';
end

% pole pitch angle
alpha_p = 2*pi/p;

alpha_v = deg2rad(alpha_v);
tmp_angle = pi/2 - alpha_v/2;

x1 = (ID/2);
y1 = 0;

x4 = (OD/2);
y4 = 0;

x5 = (OD/2) * cos(alpha_p/2);
y5 = (OD/2) * sin(alpha_p/2);

x6 = (ID/2) * cos(alpha_p/2);
y6 = (ID/2) * sin(alpha_p/2);

% finding circle center
tmp_a = OD/2 - wtrib - dm/2;
tmp_d = bm0/2 + dm/2;
tmp_m = tan(alpha_p/2);
x_sol = roots([1+tmp_m^2,-2*tmp_m*tmp_d/cos(alpha_p/2),-tmp_a^2+tmp_d^2/cos(alpha_p/2)^2]);

x9 = max(x_sol);
y9 = -tmp_d/cos(alpha_p/2) + tmp_m*x9;

[x7,y7] = g.getIntersectionRayCircle(x9,y9,-cos(tmp_angle),sin(tmp_angle),x9,y9,dm/2);
[x8,y8] = g.getIntersectionRayCircle(x9,y9,cos(tmp_angle),-sin(tmp_angle),x9,y9,dm/2);
[x3,y3] = g.getIntersectionLineLine(x8,y8,cos(alpha_v/2),sin(alpha_v/2),0,wrrib/2,1,0);
[x2,y2] = g.getIntersectionLineLine(x7,y7,cos(alpha_v/2),sin(alpha_v/2),x3,y3,cos(tmp_angle),-sin(tmp_angle));

o = g.addPoint(0,0);
p1 = g.addPoint(x1,y1);
p2 = g.addPoint(x2,y2);
p3 = g.addPoint(x3,y3);
p4 = g.addPoint(x4,y4);
p5 = g.addPoint(x5,y5);
p6 = g.addPoint(x6,y6);
p7 = g.addPoint(x7,y7);
p8 = g.addPoint(x8,y8);
p9 = g.addPoint(x9,y9);
p10 = g.addPoint(x2,y3);

if wrrib == 0

    e1 = g.addSegment(p1,p10);
    e2 = g.addSegment(p10,p3);
    e3 = g.addSegment(p3,p4);
    e4 = g.addArc(o,p4,p5,1);
    e5 = g.addSegment(p5,p6);
    e6 = g.addArc(o,p6,p1,0);
    e7 = g.addSegment(p2,p3);
    e8 = g.addSegment(p3,p8);
    e9 = g.addSegment(p8,p7);
    e10 = g.addSegment(p7,p2);
    e11 = g.addArc(p9,p8,p7,1);
    e12 = g.addSegment(p2,p10);

    l1 = g.addLoop(e1,-e12,-e10,-e11,-e8,e3,e4,e5,e6);
    l2 = g.addLoop(e7,e8,e9,e10);
    l3 = g.addLoop(-e9,e11);
    l4 = g.addLoop(e2,-e7,e12);

    g.addFace(name1, l1);
    g.addFace(name2, l2);
    g.addFace(name3 + "1", l3);
    g.addFace(name3 + "2", l4);

else

    e1 = g.addSegment(p1,p4);
    e2 = g.addArc(o,p4,p5,1);
    e3 = g.addSegment(p5,p6);
    e4 = g.addArc(o,p6,p1,0);
    e5 = g.addSegment(p2,p3);
    e6 = g.addSegment(p3,p8);
    e7 = g.addSegment(p8,p7);
    e8 = g.addSegment(p7,p2);
    e9 = g.addArc(p9,p8,p7,1);
    e10 = g.addSegment(p2,p10);
    e11 = g.addSegment(p10,p3);

    l1 = g.addLoop(e1,e2,e3,e4);
    l2 = g.addLoop(e11,e6,e9,e8,e10);
    l3 = g.addLoop(e5,e6,e7,e8);
    l4 = g.addLoop(e9,-e7);
    l5 = g.addLoop(e11,-e5,e10);

    g.addFace(name1, l1,l2);
    g.addFace(name2, l3);
    g.addFace(name3 + "1", l4);
    g.addFace(name3 + "2", l5);

end

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,160,78,146);
g.setFaceColor(name3 + "1",0,255,255);
g.setFaceColor(name3 + "2",0,255,255);

% visualizations for debug
close all;
g.setMeshMaxLength(0.2);
if nargin == 0, g.showSketch; end
if nargin == 0, g.showFaces; end

end
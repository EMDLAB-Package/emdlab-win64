% rotor & magnet
% outer rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_rm_ipm15(g, ID, OD, p, dm, wsm, w1, w2, w3, w4, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    ID = 40;
    OD = 100;
    p = 6;
    dm = 4;
    wsm = 10;
    w1 = 4;
    w2 = 0.6;
    w3 = 1;
    w4 = 1;
    name1 = 'rotor';
    name2 = 'magnet';
    name3 = 'rap';
end

% pole pitch angle
alpha_p = 2*pi/p;
tmp_angle1 = pi/2 - alpha_p/2;
tmp_angle2 = pi/2 - alpha_p;

x1 = (ID/2);
y1 = 0;

x4 = (OD/2);
y4 = 0;

x5 = (OD/2) * cos(alpha_p/2);
y5 = (OD/2) * sin(alpha_p/2);

x6 = (ID/2) * cos(alpha_p/2);
y6 = (ID/2) * sin(alpha_p/2);

% finding circle center
tmp_a = OD/2 - w4 - dm/2;
tmp_d = w1/2 + dm/2;
tmp_m = tan(alpha_p/2);
x_sol = roots([1+tmp_m^2,-2*tmp_m*tmp_d/cos(alpha_p/2),-tmp_a^2+tmp_d^2/cos(alpha_p/2)^2]);

x17 = max(x_sol);
y17 = -tmp_d/cos(alpha_p/2) + tmp_m*x17;

[x7,y7] = g.getIntersectionRayCircle(x17,y17,cos(tmp_angle1),-sin(tmp_angle1),x17,y17,dm/2);
[x8,y8] = g.getIntersectionRayCircle(x17,y17,-cos(tmp_angle1),sin(tmp_angle1),x17,y17,dm/2);

ux = cos(alpha_p/2);
uy = sin(alpha_p/2);

x9 = x8 - wsm*ux;
y9 = y8 - wsm*uy;

x10 = x7 - wsm*ux;
y10 = y7 - wsm*uy;

x14 = x10 - w2*ux;
y14 = y10 - w2*uy;

x15 = x9 - w2*ux;
y15 = y9 - w2*uy;

[ux1,uy1] = emdlab_g2d_rotatePointsXY(-ux,-uy,tmp_angle2);

x13 = x14 + w3*ux1;
y13 = y14 + w3*uy1;

x16 = x13 - dm;
y16 = y13;

x11 = x13;
y11 = y13 - w2;

x12 = x16;
y12 = y16 - w2;

x2 = x12;
y2 = 0;

x3 = x11;
y3 = 0;

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
p11 = g.addPoint(x11,y11);
p12 = g.addPoint(x12,y12);
p13 = g.addPoint(x13,y13);
p14 = g.addPoint(x14,y14);
p15 = g.addPoint(x15,y15);
p16 = g.addPoint(x16,y16);
p17 = g.addPoint(x17,y17);
o = g.addPoint(0,0);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addSegment(p3,p4);
e4 = g.addArc(o,p4,p5,1);
e5 = g.addSegment(p5,p6);
e6 = g.addArc(o,p6,p1,0);
e7 = g.addSegment(p3,p11);
e8 = g.addSegment(p11,p12);
e9 = g.addSegment(p12,p2);
e10 = g.addSegment(p13,p14);
e11 = g.addSegment(p14,p15);
e12 = g.addSegment(p15,p16);
e13 = g.addSegment(p16,p13);
e14 = g.addSegment(p7,p8);
e15 = g.addSegment(p8,p9);
e16 = g.addSegment(p9,p10);
e17 = g.addSegment(p10,p7);
e18 = g.addArc(p17,p7,p8,1);

l1 = g.addLoop(e1,-e9,-e8,-e7,e3,e4,e5,e6);
l2 = g.addLoop(e10,e11,e12,e13);
l3 = g.addLoop(e17,e18,e15,e16);
l4 = g.addLoop(e2,e7,e8,e9);
l5 = g.addLoop(e14,e15,e16,e17);
l6 = g.addLoop(e18,-e14);

g.addFace(name1, l1, l2, l3);
g.addFace(name2 + "1", l4);
g.addFace(name2 + "2", l5);
g.addFace(name3 + "1", l2);
g.addFace(name3 + "2", l6);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2 + "1",160,78,146);
g.setFaceColor(name2 + "2",160,78,146);
g.setFaceColor(name3 + "1",0,255,255);
g.setFaceColor(name3 + "2",0,255,255);

% visualizations for debug
close all;
g.setMeshMaxLength(0.2);
if nargin == 0, g.showSketch; end
if nargin == 0, g.showFaces; end

end
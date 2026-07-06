% rotor & magnet
% outer rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_rm_ipm17(g, ID, OD, p, dm1, dm2, alpha_v, g_v, wtrib, wrrib, d0, d1, d2, hD, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    ID = 90;
    OD = 150;
    p = 8;
    dm1 = 4;
    dm2 = 3;    
    alpha_v = 130;
    g_v = 1.1;
    wtrib = 1;
    wrrib = 1;
    d0 = 6;
    d1 = 5;
    d2 = 3;
    hD = 8;
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
tmp_a = OD/2 - wtrib - dm1/2;
tmp_d = d0/2 + dm1/2;
tmp_m = tan(alpha_p/2);
x_sol = roots([1+tmp_m^2,-2*tmp_m*tmp_d/cos(alpha_p/2),-tmp_a^2+tmp_d^2/cos(alpha_p/2)^2]);

x9 = max(x_sol);
y9 = -tmp_d/cos(alpha_p/2) + tmp_m*x9;

[x7,y7] = g.getIntersectionRayCircle(x9,y9,-cos(tmp_angle),sin(tmp_angle),x9,y9,dm1/2);
[x8,y8] = g.getIntersectionRayCircle(x9,y9,cos(tmp_angle),-sin(tmp_angle),x9,y9,dm1/2);
[x3,y3] = g.getIntersectionLineLine(x8,y8,cos(alpha_v/2),sin(alpha_v/2),0,wrrib/2,1,0);
[x2,y2] = g.getIntersectionLineLine(x7,y7,cos(alpha_v/2),sin(alpha_v/2),x3,y3,cos(tmp_angle),-sin(tmp_angle));

alpha_v2 = alpha_v*g_v;
alpha_v2 = min(alpha_v2,deg2rad(179));
alpha_v2 = max(alpha_v2,alpha_v);
tmp_angle2 = pi/2 - alpha_v2/2;

% finding circle center
ux = cos(alpha_v/2);
uy = sin(alpha_v/2);
m = uy/ux;
tmp_a = OD/2 - wtrib - dm2/2;
tmp_b = y3 - (0.5*dm2 + d1 + x3*uy)/ux;
x_sol = roots([1+m^2,2*m*tmp_b,tmp_b^2-tmp_a^2]);

x16 = max(x_sol);
y16 = sqrt(tmp_a^2 - x16^2);

[x14,y14] = g.getIntersectionRayCircle(x16,y16,-cos(tmp_angle2),sin(tmp_angle2),x16,y16,dm2/2);
[x13,y13] = g.getIntersectionRayCircle(x16,y16,cos(tmp_angle2),-sin(tmp_angle2),x16,y16,dm2/2);
[x12,y12] = g.getIntersectionLineLine(x13,y13,cos(alpha_v2/2),sin(alpha_v2/2),0,wrrib/2,1,0);
[x15,y15] = g.getIntersectionLineLine(x14,y14,cos(alpha_v2/2),sin(alpha_v2/2),x12,y12,cos(tmp_angle2),-sin(tmp_angle2));

x11 = x15;
y11 = y12;

x17 = (ID/2 + d2) * cos(alpha_p/2);
y17 = (ID/2 + d2) * sin(alpha_p/2);

x18 = (ID/2 + d2 + hD/2) * cos(alpha_p/2);
y18 = (ID/2 + d2 + hD/2) * sin(alpha_p/2);

x19 = (ID/2 + d2 + hD) * cos(alpha_p/2);
y19 = (ID/2 + d2 + hD) * sin(alpha_p/2);

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
p11 = g.addPoint(x11,y11);
p12 = g.addPoint(x12,y12);
p13 = g.addPoint(x13,y13);
p14 = g.addPoint(x14,y14);
p15 = g.addPoint(x15,y15);
p16 = g.addPoint(x16,y16);
p17 = g.addPoint(x17,y17);
p18 = g.addPoint(x18,y18);
p19 = g.addPoint(x19,y19);
o = g.addPoint(0,0);

if wrrib == 0

    e1 = g.addSegment(p1,p10);
    e2 = g.addSegment(p10,p3);
    e3 = g.addSegment(p3,p11);
    e4 = g.addSegment(p11,p12);
    e5 = g.addSegment(p12,p4);
    e6 = g.addArc(o,p4,p5,1);
    e7 = g.addSegment(p5,p19);
    e8 = g.addArc(o,p6,p1,0);
    e9 = g.addSegment(p2,p3);
    e10 = g.addSegment(p3,p8);
    e11 = g.addSegment(p8,p7);
    e12 = g.addSegment(p7,p2);
    e13 = g.addArc(p9,p8,p7,1);
    e14 = g.addSegment(p2,p10);
    e15 = g.addSegment(p12,p13);
    e16 = g.addSegment(p13,p14);
    e17 = g.addSegment(p14,p15);
    e18 = g.addSegment(p15,p12);
    e19 = g.addArc(p16,p13,p14,1);
    e20 = g.addSegment(p15,p11);
    e21 = g.addSegment(p19,p17);
    e22 = g.addSegment(p17,p6);
    e23 = g.addArc(p18,p17,p19,1);

    l1 = g.addLoop(e1,-e14,-e12,-e13,-e10,e3,-e20,-e17,-e19,-e15,e5,e6,e7,-e23,e22,e8);
    l2 = g.addLoop(e9,e10,e11,e12);
    l3 = g.addLoop(e15,e16,e17,e18);
    l4 = g.addLoop(e13,-e11);
    l5 = g.addLoop(e2,-e9,e14);
    l6 = g.addLoop(e19,-e16);
    l7 = g.addLoop(e4,-e18,e20);
    l8 = g.addLoop(e23,e21);

    g.addFace(name1, l1);
    g.addFace(name2 + "1", l2);
    g.addFace(name2 + "2", l3);
    g.addFace(name3 + "1", l4);
    g.addFace(name3 + "2", l5);
    g.addFace(name3 + "3", l6);
    g.addFace(name3 + "4", l7);
    g.addFace(name3 + "5", l8);

else

    e1 = g.addSegment(p1,p4);
    e2 = g.addArc(o,p4,p5,1);
    e3 = g.addSegment(p5,p19);
    e4 = g.addArc(o,p6,p1,0);
    e5 = g.addSegment(p2,p3);
    e6 = g.addSegment(p3,p8);
    e7 = g.addSegment(p8,p7);
    e8 = g.addSegment(p7,p2);
    e9 = g.addArc(p9,p8,p7,1);
    e10 = g.addSegment(p2,p10);
    e11 = g.addSegment(p10,p3);
    e12 = g.addSegment(p11,p12);
    e13 = g.addSegment(p12,p13);
    e14 = g.addSegment(p13,p14);
    e15 = g.addSegment(p14,p15);
    e16 = g.addSegment(p15,p11);
    e17 = g.addSegment(p15,p12);
    e18 = g.addArc(p16,p13,p14,1);
    e19 = g.addSegment(p19,p17);
    e20 = g.addSegment(p17,p6);
    e21 = g.addArc(p18,p17,p19,1);

    e22 = g.splitArc(e9);
    e23 = g.splitArc(e18);
    e24 = g.splitArc(e21);

    l1 = g.addLoop(e1,e2,e3,-e24,-e21,e20,e4);
    l2 = g.addLoop(e11,e6,e9,e22,e8,e10);
    l3 = g.addLoop(e12,e13,e18,e23,e15,e16);
    l4 = g.addLoop(e5,e6,e7,e8);
    l5 = g.addLoop(e13,e14,e15,e17);
    l6 = g.addLoop(e9,e22,-e7);
    l7 = g.addLoop(e11,-e5,e10);
    l8 = g.addLoop(e18,e23,-e14);
    l9 = g.addLoop(e12,-e17,e16);
    l10 = g.addLoop(e21,e24,e19);

    g.addFace(name1, l1, l2, l3);
    g.addFace(name2 + "1", l4);
    g.addFace(name2 + "2", l5);
    g.addFace(name3 + "1", l6);
    g.addFace(name3 + "2", l7);
    g.addFace(name3 + "3", l8);
    g.addFace(name3 + "4", l9);
    g.addFace(name3 + "5", l10);

end

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2 + "1",160,78,146);
g.setFaceColor(name2 + "2",160,78,146);
g.setFaceColor(name3 + "1",0,255,255);
g.setFaceColor(name3 + "2",0,255,255);
g.setFaceColor(name3 + "3",0,255,255);
g.setFaceColor(name3 + "4",0,255,255);
g.setFaceColor(name3 + "5",0,255,255);

% visualizations for debug
close all;
g.setMeshMaxLength(0.2);
if nargin == 0, g.showSketch; end
if nargin == 0, g.showFaces; end

end
% rotor & magnet
% inner rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_tc_wfsm1(g, Dsh, ISD, gap, gg, Nagl, p, embrace, wpole, hs0, ghc, gdc, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    Dsh = 56;
    ISD = 150;
    gap = 1;
    gg = 3;
    Nagl = 2;
    p = 6;
    embrace = 0.73;
    wpole = 25;
    hs0 = 5;
    ghc = 0.8;
    gdc = 0.8;
    name1 = 'rotor';
    name2 = 'fcoil';
    name3 = 'rap';
end

% rotor pole pitch angle
alpha_p = 2*pi/p;

% limit embrace
embrace = min(embrace,0.95);
embrace = max(embrace,0.05);

% magnet arc angle
alpha_t = embrace * alpha_p;

x1 = Dsh/2;
y1 = 0;

x2 = ISD/2 - gap;
y2 = 0;

x3 = x2 + gap/Nagl;
y3 = 0;

x13 = (ISD/2 - gg*gap) * cos(alpha_t/2);
y13 = (ISD/2 - gg*gap) * sin(alpha_t/2);

u = [x13,y13];
u = u/norm(u);

x12 = x13 - hs0*u(1);
y12 = y13 - hs0*u(2);

x9 = x12;
y9 = wpole/2;

x10 = x9;
y10 = y9 + ghc * (y12 - y9);

[xi,~] = g.getIntersectionLineLine(x10,y10,1,0,0,0,cos(alpha_p/2),sin(alpha_p/2));

x8 = x9 - gdc * (x9-xi);
y8 = y9;

x11 = x8;
y11 = y10;

x4 = x3 * cos(alpha_p/2);
y4 = x3 * sin(alpha_p/2);

r_tmp = norm([x13,y13]);
x5 = r_tmp * cos(alpha_p/2);
y5 = r_tmp * sin(alpha_p/2);

r_tmp = norm([x11,y11]);
x6 = r_tmp * cos(alpha_p/2);
y6 = r_tmp * sin(alpha_p/2);

x7 = x1 * cos(alpha_p/2);
y7 = x1 * sin(alpha_p/2);

xc = (x13^2+y13^2-x2^2)/(2*(x13-x2));

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
o = g.addPoint(0,0);
o2 = g.addPoint(xc,0);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addArc(o,p3,p4,1);
e4 = g.addSegment(p4,p5);
e5 = g.addSegment(p5,p6);
e6 = g.addSegment(p6,p7);
e7 = g.addArc(o,p7,p1,0);
e8 = g.addArc(o2,p2,p13,1);
e9 = g.addSegment(p13,p12);
e10 = g.addSegment(p12,p10);
e11 = g.addSegment(p8,p9);
e12 = g.addSegment(p9,p10);
e13 = g.addSegment(p10,p11);
e14 = g.addSegment(p11,p8);
e15 = g.addArc(o,p11,p6,1);
e16 = g.addArc(o,p5,p13,0);

% add loops
l1 = g.addLoop(e1,e8,e9,e10,-e12,-e11,-e14,e15,e6,e7);
l2 = g.addLoop(e11,e12,e13,e14);
l3 = g.addLoop(e2,e3,e4,e16,-e8);
l4 = g.addLoop(e5,-e15,-e13,-e10,-e9,-e16);

% add faces
g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3 + "1", l3);
g.addFace(name3 + "2", l4);

% set face colors
g.setFaceColor(name1,200,200,200)
g.setFaceColor(name2,160,78,146);
g.setFaceColor(name3 + "1",0,255,255)
g.setFaceColor(name3 + "2",0,255,255)

% visualizations for debug
close all;
g.setMeshMaxLength(0.2);
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
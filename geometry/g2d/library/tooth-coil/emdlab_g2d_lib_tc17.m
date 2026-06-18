% EMDLAB: Electrical Machines Design Laboratory
% tooth & coil geometry template

function emdlab_g2d_lib_tc17(g, ISD, OSD, Ns, wst, dss, bs0, hs0, tta, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    ISD = 92;
    OSD = 140;
    Ns = 12;
    wst = 13;
    dss = 16;
    bs0 = 5;
    hs0 = 2;
    tta = 10;
    name1 = 'stator';
    name2 = 'sca';
    name3 = 'sap';
end

% check unfeasible geometries
if (ISD/2+dss) > OSD/2
    error('OD/2 must be higher than (ID/2+ds)');
end

% stator slot pitch angle
alpha_s = 2*pi/Ns;

% slot openning angle
gamma_so = 2*asin(bs0/ISD);

% tooth tip angle in radian
tta =deg2rad(tta);

x1 = (ISD/2) * cos(gamma_so/2);
y1 = (ISD/2) * sin(gamma_so/2);

x2 = (ISD/2+hs0) * cos(gamma_so/2);
y2 = (ISD/2+hs0) * sin(gamma_so/2);

u = [x2,y2];
u = emdlab_g2d_rotatePoints(u,pi/2-tta);

[x3,y3] = g.getIntersectionLineLine(x2,y2,u(1),u(2),0,-wst*0.5/cos(alpha_s/2),cos(alpha_s/2),sin(alpha_s/2));

x5 = ISD/2 + dss;
y5 = 0;

[x4,y4] = g.getIntersectionLineLine(x5,y5,-sin(alpha_s/2),cos(alpha_s/2),0,-wst*0.5/cos(alpha_s/2),cos(alpha_s/2),sin(alpha_s/2));

x6 = OSD/2;
y6 = 0;

x7 = (OSD/2) * cos(alpha_s/2);
y7 = (OSD/2) * sin(alpha_s/2);

x8 = (ISD/2) * cos(alpha_s/2);
y8 = (ISD/2) * sin(alpha_s/2);

x9 = (ISD/2);
y9 = 0;

[x10,y10] = g.getIntersectionLineLine(x3,y3,-sin(alpha_s/2),cos(alpha_s/2),0,0,1,0);

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
p10 = g.addPoint(x10,y10);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addSegment(p3,p4);
e4 = g.addSegment(p4,p5);
e5 = g.addSegment(p5,p6);
e6 = g.addArc(o,p6,p7,1);
e7 = g.addSegment(p7,p8);
e8 = g.addArc(o,p8,p1,0);
e9 = g.addSegment(p9,p10);
e10 = g.addSegment(p10,p5);
e11 = g.addArc(o,p9,p1,1);
e12 = g.addSegment(p10,p3);

l1 = g.addLoop(e1,e2,e3,e4,e5,e6,e7,e8);
l2 = g.addLoop(e10,-e4,-e3,-e12);
l3 = g.addLoop(e9,e12,-e2,-e1,-e11);

g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3, l3);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,255,137,39);
g.setFaceColor(name3,0,255,255);

% visualizations for debug
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
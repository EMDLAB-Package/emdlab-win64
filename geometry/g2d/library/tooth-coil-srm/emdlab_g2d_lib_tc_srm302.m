function emdlab_g2d_lib_tc_srm302(g, ID, OD, Ntooth, g_beta, g_wy, name1, name2)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    ID = 45;
    OD = 74;
    Ntooth = 6;
    g_beta = 0.4;
    g_wy = 0.4;
    name1 = 'rotor';
    name2 = 'rca';
end

% calculate dependent variables.
alpha_t = 2*pi/Ntooth;
beta_t = g_beta * alpha_t;
wt = 2 * (OD/2) * sin(beta_t/2);
wy = g_wy * wt;
gamma_tr = 2*asin(wt*0.5/(ID/2+wy));

% adding points
x1 = ID/2;
y1 = 0;

x2 = ID/2+wy;
y2 = 0;

x3 = OD/2;
y3 = 0;

x4 = (OD/2) * cos(beta_t/2);
y4 = (OD/2) * sin(beta_t/2);

x5 = (ID/2+wy) * cos(gamma_tr/2);
y5 = (ID/2+wy) * sin(gamma_tr/2);

x6 = (ID/2+wy) * cos(alpha_t/2);
y6 = (ID/2+wy) * sin(alpha_t/2);

x7 = (ID/2) * cos(alpha_t/2);
y7 = (ID/2) * sin(alpha_t/2);

x8 = (ID/2) * cos(gamma_tr/2);
y8 = (ID/2) * sin(gamma_tr/2);

x9 = (OD/2) * cos(alpha_t/2);
y9 = (OD/2) * sin(alpha_t/2);

% adding points
p1 = g.addPoint(x1,y1);
p2 = g.addPoint(x2,y2);
p3 = g.addPoint(x3,y3);
p4 = g.addPoint(x4,y4);
p5 = g.addPoint(x5,y5);
p6 = g.addPoint(x6,y6);
p7 = g.addPoint(x7,y7);
p8 = g.addPoint(x8,y8);
p9 = g.addPoint(x9,y9);
o = g.addPoint(0,0);

% adding edges
e1 = g.addSegment(p1, p2);
e2 = g.addSegment(p2, p3);
e3 = g.addArc(o, p3, p4, 1);
e4 = g.addSegment(p4, p5);
e5 = g.addArc(o, p5, p6, 1);
e6 = g.addSegment(p6, p7);
e7 = g.addArc(o, p7, p8, 0);
e8 = g.addArc(o, p8, p1, 0);
e9 = g.addArc(o, p2, p5, 1);
e10 = g.addSegment(p5, p8);
e11 = g.addSegment(p9, p6);
e12 = g.addArc(o, p4, p9, 1);

% adding loops
l1 = g.addLoop(e1,e9,e10,e8);
l2 = g.addLoop(e2,e3,e4,-e9);
l3 = g.addLoop(e5,e6,e7,-e10);
l4 = g.addLoop(e12,e11,-e5,-e4);

% adding faces
g.addFace(name1 + "1", l1);
g.addFace(name1 + "2", l2);
g.addFace(name1 + "3", l3);
g.addFace(name2, l4);

% set default colors
g.setFaceColor(name1 + "1",200,200,200);
g.setFaceColor(name1 + "2",200,200,200);
g.setFaceColor(name1 + "3",200,200,200);
g.setFaceColor(name2,255,137,39);

% visualizations for debug
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
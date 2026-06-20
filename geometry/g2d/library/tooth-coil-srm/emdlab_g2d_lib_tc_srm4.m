% with rectangular coil

function emdlab_g2d_lib_tc_srm4(g, ID, OD, Ntooth, g_beta, g_wy, name1, name2, name3)

if nargin == 0
    g = emdlab_g2d_db;
    ID = 75;
    OD = 125;
    Ntooth = 8;
    g_beta = 0.3;
    g_wy = 0.4;
    name1 = 'rotor';
    name2 = 'rc';
    name3 = 'rap';
end

% calculate dependent variables.
alpha_t = 2*pi/Ntooth;
beta_t = g_beta * alpha_t;
wt = 2 * (OD/2) * sin(beta_t/2);
wy = g_wy * wt;
gamma_tr = 2*asin(wt*0.5/(ID/2+wy));

% adding points
p1 = g.addPoint(0,0);
p2 = g.addPoint(ID/2,0);
p3 = g.addPoint(OD/2,0);
p4 = g.addPoint((OD/2)*cos(beta_t/2),(OD/2)*sin(beta_t/2));
p5 = g.addPoint((OD/2)*cos(alpha_t/2),(OD/2)*sin(alpha_t/2));
[p6,p6h] = g.addPoint((ID/2+wy)*cos(alpha_t/2),(ID/2+wy)*sin(alpha_t/2));
p7 = g.addPoint((ID/2)*cos(alpha_t/2),(ID/2)*sin(alpha_t/2));
[p8,p8h] = g.addPoint((ID/2+wy)*cos(gamma_tr/2),(ID/2+wy)*sin(gamma_tr/2));

[xi,yi] = g.getIntersectionLineLine(p8h.x,p8h.y,0,1,0,0,p6h.x,p6h.y);
p9 = g.addPoint(xi,yi);
[xi,yi] = g.getIntersectionLineArc(xi,yi,1,0,0,0,OD/2,beta_t*90/pi,alpha_t*90/pi);
p10 = g.addPoint(xi,yi);
[xi,yi] = g.getIntersectionLineLine(p8h.x,p8h.y,1,0,xi,yi,0,1);
p11 = g.addPoint(xi,yi);

% adding edges
e1 = g.addSegment(p2, p3);
e2 = g.addArc(p1, p3, p4, 1);
e3 = g.addSegment(p4, p11);
e4 = g.addSegment(p11, p8);
e5 = g.addArc(p1, p8, p6, 1);
e6 = g.addSegment(p6, p7);
e7 = g.addArc(p1, p7, p2, 0);
e8 = g.addArc(p1, p4, p10, 1);
e9 = g.addArc(p1, p10, p5, 1);
e10 = g.addSegment(p5, p9);
e11 = g.addSegment(p9, p6);
e12 = g.addSegment(p11, p10);
e13 = g.addSegment(p10, p9);
e14 = g.addSegment(p9, p8);

% adding loops
l1 = g.addLoop(e1,e2,e3,e4,e5,e6,e7);
l2 = g.addLoop(e12,e13,e14,-e4);
l3 = g.addLoop(e8,-e12,-e3);
l4 = g.addLoop(e9,e10,-e13);
l5 = g.addLoop(e11,-e5,-e14);

% adding faces
g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3+"1", l3);
g.addFace(name3+"2", l4);
g.addFace(name3+"3", l5);

% set default colors
g.setFaceColor(name1,200,200,200)
g.setFaceColor(name2,255,137,39);
g.setFaceColor(name3+"1",0,255,255);
g.setFaceColor(name3+"2",0,255,255);
g.setFaceColor(name3+"3",0,255,255);

% visualize for debug
if nargin ==0, g.showSketch; end

end
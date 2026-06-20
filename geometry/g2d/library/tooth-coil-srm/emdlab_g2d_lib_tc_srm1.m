function emdlab_g2d_lib_tc_srm1(g, ID, OD, Ntooth, g_beta, g_wy, name1, name2)

% calculate dependent variables.
alpha_t = 2*pi/Ntooth;
beta_t = g_beta * alpha_t;
wt = 2 * (ID/2) * sin(beta_t/2);
wy = g_wy * wt;
gamma_tr = 2*asin(wt*0.5/(OD/2-wy));

% adding points
p1 = g.addPoint(0,0);
p2 = g.addPoint(ID/2,0);
p3 = g.addPoint(OD/2,0);
p4 = g.addPoint((OD/2)*cos(alpha_t/2),(OD/2)*sin(alpha_t/2));
p5 = g.addPoint((OD/2-wy)*cos(alpha_t/2),(OD/2-wy)*sin(alpha_t/2));
p6 = g.addPoint((ID/2)*cos(alpha_t/2),(ID/2)*sin(alpha_t/2));
p7 = g.addPoint((ID/2)*cos(beta_t/2),(ID/2)*sin(beta_t/2));
p8 = g.addPoint((OD/2-wy)*cos(gamma_tr/2),(OD/2-wy)*sin(gamma_tr/2));

% adding edges
e1 = g.addSegment(p2, p3);
e2 = g.addArc(p1, p3, p4, 1);
e3 = g.addSegment(p4, p5);
e4 = g.addSegment(p5, p6);
e5 = g.addArc(p1, p6, p7, 0);
e6 = g.addArc(p1, p7, p2, 0);
e7 = g.addSegment(p7, p8);
e8 = g.addArc(p1, p8, p5, 1);

% adding loops
l1 = g.addLoop(e1,e2,e3,-e8,-e7,e6);
l2 = g.addLoop(e7,e8,e4,e5);

% adding faces
g.addFace(name1, l1);
g.addFace(name2, l2);

% set default colors
g.setFaceColor(name1,200,200,200)
g.setFaceColor(name2,255,137,39);

end
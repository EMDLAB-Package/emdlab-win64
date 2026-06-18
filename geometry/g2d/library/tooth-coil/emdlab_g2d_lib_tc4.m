% stator tooth & coil

function emdlab_g2d_lib_tc4(g, ID, OD, Ns, wst, dss, bs0, hs0, tta, name1, name2, name3)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    ID (1,1) double {mustBePositive}
    OD (1,1) double {mustBePositive}
    Ns (1,1) double {mustBePositive,mustBeInteger}
    wst (1,1) double {mustBePositive}
    dss (1,1) double {mustBePositive}
    bs0 (1,1) double {mustBePositive}
    hs0 (1,1) double {mustBePositive}
    tta (1,1) double {mustBePositive}
    name1 (1,:) char;
    name2 (1,:) char;
    name3 (1,:) char;

end

gamma_so = 2*asin(bs0/(OD-2*hs0));
alpha_s = 2*pi/Ns;
tta = tta * pi/180;

x1 = (OD/2) * cos(gamma_so/2);
y1 = (OD/2) * sin(gamma_so/2);

u = [x1,y1];
u = u/norm(u);

x2 = x1 - hs0 * u(1);
y2 = y1 - hs0 * u(2);

[ux,uy] = emdlab_g2d_rotatePointsXY(u(1), u(2), pi/2+tta);
u1 = [ux,uy];

u2 = [cos(alpha_s/2), sin(alpha_s/2)];

[x3,y3] = emdlab_g2d_getLineLineIntersection([0,-wst*0.5/cos(alpha_s/2)], u2, [x2,y2], u1);

[ux,uy] = emdlab_g2d_rotatePointsXY(u2(1), u2(2), -pi/2);
u3 = [ux,uy];

    function err = err_fcn(r)
        [xtmp,ytmp] = emdlab_g2d_getLineLineIntersection([0,-wst*0.5/cos(alpha_s/2)], u2, [OD/2-dss+r,0], u3);
        r2 = norm([xtmp,ytmp]-[OD/2-dss+r,0]);
        err = r2-r;
    end

r = fzero(@err_fcn,1);
[x4,y4] = emdlab_g2d_getLineLineIntersection([0,-wst*0.5/cos(alpha_s/2)], u2, [OD/2-dss+r,0], u3);

x5 = OD/2-dss;
y5 = 0;

x6 = ID/2;
y6 = 0;

[x7,y7] = emdlab_g2d_rotatePointsXY(x6,y6,alpha_s/2);

[x8,y8] = emdlab_g2d_rotatePointsXY(x1,y1,alpha_s/2-gamma_so/2);

x9 = OD/2;
y9 = 0;

x10 = sqrt(x3^2+y3^2);
y10 = 0;

x11 = OD/2-dss+r;
y11 = 0;

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
p11 = g.addPoint(x11,y11);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addSegment(p3,p4);
e4 = g.addArc(p11,p4,p5,1);
e5 = g.addSegment(p5,p6);
e6 = g.addArc(o,p6,p7,1);
e7 = g.addSegment(p7,p8);
e8 = g.addArc(o,p8,p1,0);
e9 = g.addSegment(p9,p10);
e10 = g.addSegment(p10,p5);
e11 = g.addArc(o,p9,p1,1);
e12 = g.addArc(o,p10,p3,1);

l1 = g.addLoop(-e8,-e7,-e6,-e5,-e4,-e3,-e2,-e1);
l2 = g.addLoop(e12,e3,e4,-e10);
l3 = g.addLoop(e11,e1,e2,-e12,-e9);

g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3, l3);

g.setFaceColor(name1,200,200,200)
g.setFaceColor(name2,255,137,39);
g.setFaceColor(name3,0,255,255)

end
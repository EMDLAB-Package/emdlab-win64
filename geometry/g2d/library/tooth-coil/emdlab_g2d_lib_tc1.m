% EMDLAB: Electrical Machines Design Laboratory
% stator tooth & coil

function emdlab_g2d_lib_tc1(g, ISD, OSD, Ns, wst, dss, bs0, hs0, tta, name1, name2, name3)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    ISD (1,1) double {mustBePositive}
    OSD (1,1) double {mustBePositive}
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

% check unfeasible geometries
if (ISD/2+dss) > OSD/2
    error('OD/2 must be higher than (ID/2+ds)');
end

gamma_so = 2*asin(bs0/ISD);
alpha_s = 2*pi/Ns;
tta = tta * pi/180;

o = emdlab_g2d_point(0,0);
p1 = emdlab_g2d_point((ISD/2) * cos(gamma_so/2), (ISD/2) * sin(gamma_so/2));
p2 = p1.getMoveRadial(hs0);
l = emdlab_g2d_line(emdlab_g2d_point(0,0),emdlab_g2d_point(cos(alpha_s/2),sin(alpha_s/2)));
l.moveNormal(-wst/2);
p3 = p2.getMoveRadial(1);
p3.rotateAroundPoint(p2,pi/2-tta);
u23 = p3 - p2;
u23.normalize;

    function d = tmp_f(x)
        p3.x = p2.x + u23.x * x;
        p3.y = p2.y + u23.y * x;
        d = p3.getDistanceFromLine(l);
    end

x = fminsearch(@tmp_f,1);
tmp_f(x);

p5 = emdlab_g2d_point(ISD/2 + dss, 0);
alpha_tmp = alpha_s/2 - asin(wst/(ISD+2*dss));

p4 = p5.getRotateAroundPoint(o,alpha_tmp);
p6 = emdlab_g2d_point(OSD/2,0);
p7 = p6.getRotateAroundPoint(o,alpha_s/2);
p8 = p7.getMoveRadial(-OSD/2+ISD/2);

p9 = emdlab_g2d_point(ISD/2,0);
tmp = p3.getDistanceFromOrigin;
p10 = emdlab_g2d_point(tmp,0);

o = g.addPoint(0,0);
p1 = g.addPoint(p1);
p2 = g.addPoint(p2);
p3 = g.addPoint(p3);
p4 = g.addPoint(p4);
p5 = g.addPoint(p5);
p6 = g.addPoint(p6);
p7 = g.addPoint(p7);
p8 = g.addPoint(p8);
p9 = g.addPoint(p9);
p10 = g.addPoint(p10);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addSegment(p3,p4);
e4 = g.addArc(o,p4,p5,0);
e5 = g.addSegment(p5,p6);
e6 = g.addArc(o,p6,p7,1);
e7 = g.addSegment(p7,p8);
e8 = g.addArc(o,p8,p1,0);
e9 = g.addSegment(p9,p10);
e10 = g.addSegment(p10,p5);
e11 = g.addArc(o,p9,p1,1);
e12 = g.addArc(o,p10,p3,1);

l1 = g.addLoop(e1,e2,e3,e4,e5,e6,e7,e8);
l2 = g.addLoop(e10,-e4,-e3,-e12);
l3 = g.addLoop(e9,e12,-e2,-e1,-e11);

g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3, l3);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,255,137,39);
g.setFaceColor(name3,0,255,255);

end
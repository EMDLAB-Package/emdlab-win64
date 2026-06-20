% rotor & magnet
% outer rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_rm_spm2(g, ID, OD, p, dm, embrace, name1, name2, name3)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    ID (1,1) double {mustBePositive}
    OD (1,1) double {mustBePositive}
    p (1,1) double {mustBePositive,mustBeInteger}
    dm (1,1) double {mustBePositive}
    embrace (1,1) double {mustBePositive}
    name1 (1,:) char;
    name2 (1,:) char;
    name3 (1,:) char;

end

alpha_p = 2*pi/p;

o = emdlab_g2d_point(0,0);
p1 = emdlab_g2d_point(ID/2,0);
p2 = emdlab_g2d_point(ID/2+dm,0);
p3 = emdlab_g2d_point(OD/2,0);
p4 = p1.getRotateAroundOrigin(embrace*alpha_p/2);
p5 = p2.getRotateAroundOrigin(embrace*alpha_p/2);
p6 = p1.getRotateAroundOrigin(alpha_p/2);
p7 = p2.getRotateAroundOrigin(alpha_p/2);
p8 = p3.getRotateAroundOrigin(alpha_p/2);

o = g.addPoint(o);
p1 = g.addPoint(p1);
p2 = g.addPoint(p2);
p3 = g.addPoint(p3);
p4 = g.addPoint(p4);
p5 = g.addPoint(p5);
p6 = g.addPoint(p6);
p7 = g.addPoint(p7);
p8 = g.addPoint(p8);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addSegment(p4,p5);
e4 = g.addSegment(p6,p7);
e5 = g.addSegment(p7,p8);

e6 = g.addArc(o,p1,p4,1);
e7 = g.addArc(o,p2,p5,1);
e8 = g.addArc(o,p4,p6,1);
e9 = g.addArc(o,p5,p7,1);
e10 = g.addArc(o,p3,p8,1);

l1 = g.addLoop(e2,e10,-e5,-e9,-e7);
l2 = g.addLoop(e1,e7,-e3,-e6);
l3 = g.addLoop(e3,e9,-e4,-e8);

g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3, l3);

% set mesh zone colors
g.setFaceColor(name1,200,200,200)
g.setFaceColor(name2,28,255,28)
g.setFaceColor(name3,0,255,255)

end
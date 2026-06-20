% rotor & magnet
% outer rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_rm_ipm8(g, ID, OD, poles, alpha, beta, dm, gamma, wtrib, wrrib, name1, name2, name3)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    ID (1,1) double {mustBePositive}
    OD (1,1) double {mustBePositive}
    poles (1,1) double {mustBePositive,mustBeInteger}
    alpha (1,1) double {mustBePositive}
    beta (1,:) double {mustBePositive}
    dm (1,:) double {mustBePositive}
    gamma (1,:) double {mustBePositive}
    wtrib (1,:) double {mustBePositive}
    wrrib (1,1) double {mustBePositive}
    name1 (1,:) char;
    name2 (1,:) char;
    name3 (1,:) char;

end

beta = beta * pi/180;

% pole pitch angle
alpha_p = 2*pi/poles;

% origin point
o = g.addPoint(0,0);

tmp = alpha*alpha_p/2;
x1 = (OD/2 - wtrib) * cos(tmp);
y1 = (OD/2 - wtrib) * sin(tmp);
m = tan(pi/2 - beta);

y2 = wrrib/2;
x2 = (y2-y1)/m + x1;

y3 = y2;
x3 = x2 + dm/cos(beta);

[x4,y4] = g.getIntersectionRayCircle(x3,y3,1,m,0,0,OD/2-wtrib);

e1 = g.addSegmentByCoordinates(x2,y2,x3,y3);

tmp = g.addSegmentByCoordinates(x3,y3,x4,y4);
tmp = g.splitEdge(tmp,1-gamma);
e2 = tmp(1);
e3 = tmp(2);

[e4,e4h] = g.extendSegmentBySegment(e2,pi/2,dm,1);
[e5,e5h] = g.extendSegmentBySegment(e3,pi/2,dm,1);

e6 = g.addArcByCoordinates(0,0,x4,y4,x1,y1,1);

e7 = g.addSegmentByCoordinates(x1,y1,e5h.ptr.p1.x,e5h.ptr.p1.y);
e8 = g.addSegmentByCoordinates(e5h.ptr.p1.x,e5h.ptr.p1.y,e4h.ptr.p1.x,e4h.ptr.p1.y);
e9 = g.addSegmentByCoordinates(e4h.ptr.p1.x,e4h.ptr.p1.y,x2,y2);

% outer rotor loop
lro = g.addAnnularSectorLoop(ID/2, OD/2, 0, alpha_p/2, 0, 0);
lri = g.addLoop(e1,e2,e3,e6,e7,e8,e9);

lm = g.addLoop(e3,e5,e8,-e4);
la1 = g.addLoop(e6,e7,-e5);
la2 = g.addLoop(e1,e2,e4,e9);

g.addFace(name1,lro,lri);
g.addFace(name2,lm);
g.addFace(name3+"1",la1);
g.addFace(name3+"2",la2);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,160,78,146);
g.setFaceColor(name3+"1",0,255,255);
g.setFaceColor(name3+"2",0,255,255);

end
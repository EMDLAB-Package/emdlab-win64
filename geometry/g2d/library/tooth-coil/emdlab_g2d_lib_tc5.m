% EMDLAB: Electrical Machines Design Laboratory
% tooth & coil geometry template: constant slot width

function emdlab_g2d_lib_tc5(g, ID, OD, Nslots, wslot, dslot, bs0, hs0, tta, name1, name2, name3)

% % input arguments check
% arguments
%     g (1,1) emdlab_g2d_db
%     ID (1,1) double {mustBePositive}
%     OD (1,1) double {mustBePositive}
%     Nslots (1,1) double {mustBePositive,mustBeInteger}
%     wslot (1,1) double {mustBePositive}
%     dslot (1,1) double {mustBePositive}
%     bs0 (1,1) double {mustBePositive}
%     hs0 (1,1) double {mustBePositive}
%     tta (1,1) double {mustBePositive}
%     name1 (1,:) char;
%     name2 (1,:) char;
%     name3 (1,:) char;
% end

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    ID = 75;
    OD = 125;
    Nslots = 36;
    wslot = 4;
    dslot = 12;
    bs0 = 2.4;
    hs0 = 0.6;
    tta = 30;
    name1 = 'stator';
    name2 = 'sca';
    name3 = 'sap';
end

% check unfeasible geometries
if (ID/2+dslot) > OD/2
    error('OD/2 must be higher than (ID/2+ds)');
end

gamma_so = 2*asin(bs0/ID);
alpha_s = 2*pi/Nslots;
tta = tta * pi/180;

p1 = g.addPoint(0,0);
[p2,p2h] = g.addPoint((ID/2) * cos(gamma_so/2), (ID/2) * sin(gamma_so/2));
[p3, p3h] = g.addPoint(p2h.x + hs0, p2h.y);
[xi,yi] = g.getIntersectionLineLine(p3h.x, p3h.y, cos(pi/2-tta), sin(pi/2-tta),0,wslot/2,1,0);
[p4,p4h] = g.addPoint(xi,yi);
p5 = g.addPoint(ID/2 + dslot, wslot/2);
p6 = g.addPoint(ID/2 + dslot, 0);
p7 = g.addPoint(OD/2,0);
p8 = g.addPoint((OD/2)*cos(alpha_s/2), (OD/2)*sin(alpha_s/2));
p9 = g.addPoint((ID/2)*cos(alpha_s/2), (ID/2)*sin(alpha_s/2));
p10 = g.addPoint(ID/2,0);
p11 = g.addPoint(p4h.x,0);

e1 = g.addSegment(p2,p3);
e2 = g.addSegment(p3,p4);
e3 = g.addSegment(p4,p5);
e4 = g.addSegment(p5,p6);
e5 = g.addSegment(p6,p7);
e6 = g.addArc(p1,p7,p8,1);
e7 = g.addSegment(p8,p9);
e8 = g.addArc(p1,p9,p2,0);
e9 = g.addArc(p1,p2,p10,0);
e10 = g.addSegment(p10,p11);
e11 = g.addSegment(p11,p6);
e12 = g.addSegment(p4,p11);

l1 = g.addLoop(e1,e2,e3,e4,e5,e6,e7,e8);
l2 = g.addLoop(e11,-e4,-e3,e12);
l3 = g.addLoop(e10,-e12,-e2,-e1,e9);

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
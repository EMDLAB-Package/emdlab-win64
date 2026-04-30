% EMDLAB: Electrical Machines Design Laboratory
% tooth & coil geometry template: constant slot width

function emdlab_g2d_lib_tc12(g, ID, OD, Nslots, dslot, ratio, name1, name2)

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
    ID = 80;
    OD = 120;
    Nslots = 18;
    dslot = 12;
    ratio = 0.5;
    name1 = 'stator';
    name2 = 'sca';
end

% check unfeasible geometries
if (OD/2-dslot) < ID/2
    error('ID/2 must be lower than (OD/2-dslot)');
end

alpha_s = 2*pi/Nslots;
gamma_so = ratio * alpha_s;

p1 = g.addPoint(0,0);
p2 = g.addPoint(OD/2,0);
p3 = g.addPoint(OD/2 - dslot, 0);
p4 = g.addPoint(ID/2, 0);
p5 = g.addPoint((ID/2)*cos(alpha_s/2), (ID/2)*sin(alpha_s/2));
p6 = g.addPoint((OD/2)*cos(alpha_s/2), (OD/2)*sin(alpha_s/2));
p7 = g.addPoint((OD/2)*cos(gamma_so/2), (OD/2)*sin(gamma_so/2));
p8 = g.addPoint((OD/2 - dslot)*cos(gamma_so/2), (OD/2 - dslot)*sin(gamma_so/2));

e1 = g.addSegment(p2,p3);
e2 = g.addSegment(p3,p4);
e3 = g.addArc(p1,p4,p5,1);
e4 = g.addSegment(p5,p6);
e5 = g.addArc(p1,p6,p7,0);
e6 = g.addSegment(p7,p8);
e7 = g.addArc(p1,p8,p3,0);
e8 = g.addArc(p1,p7,p2,0);

l1 = g.addLoop(-e7,-e6,-e5,-e4,-e3,-e2);
l2 = g.addLoop(-e1,-e8,e6,e7);

g.addFace(name1, l1);
g.addFace(name2, l2);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,255,137,39);

% visualizations for debug
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
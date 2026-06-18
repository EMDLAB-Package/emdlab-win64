% EMDLAB: Electrical Machines Design Laboratory
% tooth & coil geometry template: constant slot width

function emdlab_g2d_lib_tc13(g, ID, OD, Nslots, dslot, wslot, name1, name2)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    ID = 80;
    OD = 120;
    Nslots = 18;
    dslot = 13;
    wslot = 7;
    name1 = 'stator';
    name2 = 'sca';
end

% check unfeasible geometries
if (ID/2+dslot) > OD/2
    error('OD/2 must be higher than (ID/2+ds)');
end

% stator slot pitch angle
alpha_s = 2*pi/Nslots;

% stator slot angle at tip
gamma_sst = 2*asin(wslot/ID);

% stator slot angle at root
gamma_ssr = 2*asin(wslot/(ID+2*dslot));

p1 = g.addPoint(0,0);
p2 = g.addPoint(ID/2,0);
p3 = g.addPoint(ID/2 + dslot, 0);
p4 = g.addPoint(OD/2, 0);
p5 = g.addPoint((OD/2)*cos(alpha_s/2), (OD/2)*sin(alpha_s/2));
p6 = g.addPoint((ID/2)*cos(alpha_s/2), (ID/2)*sin(alpha_s/2));
p7 = g.addPoint((ID/2)*cos(gamma_sst/2), (ID/2)*sin(gamma_sst/2));
p8 = g.addPoint((ID/2 + dslot)*cos(gamma_ssr/2), (ID/2 + dslot)*sin(gamma_ssr/2));

e1 = g.addSegment(p2,p3);
e2 = g.addSegment(p3,p4);
e3 = g.addArc(p1,p4,p5,1);
e4 = g.addSegment(p5,p6);
e5 = g.addArc(p1,p6,p7,0);
e6 = g.addSegment(p7,p8);
e7 = g.addArc(p1,p8,p3,0);
e8 = g.addArc(p1,p7,p2,0);

l1 = g.addLoop(e2,e3,e4,e5,e6,e7);
l2 = g.addLoop(e1,-e7,-e6,e8);

g.addFace(name1, l1);
g.addFace(name2, l2);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,255,137,39);

% visualizations for debug
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
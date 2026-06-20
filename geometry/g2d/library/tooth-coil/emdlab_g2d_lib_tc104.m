% stator tooth & coil

function emdlab_g2d_lib_tc104(g, yShift, wst, wss, dss, dsy, bs0, hs0, tta, paperT, bw0, name1, name2, name3)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    yShift = 1;
    wst = 10;
    wss = 18;
    dss = 20;
    dsy = 10;
    bs0 = 3;
    hs0 = 2;
    tta = 10;
    paperT = 1;
    bw0 = 2;
    name1 = 'stator';
    name2 = 'sca';
    name3 = 'sap';
end

tta = tta*pi/180;
ltmp = wss/2-bs0/2;
hs1 = ltmp*tan(tta);
bw0 = max(bw0,paperT);

s1 = g.addSegmentByCoordinates(0,yShift,wss/2+wst/2-bs0/2,yShift);
s2 = g.extendSegmentBySegment(s1,pi/2,hs0);
s3 = g.extendSegmentBySegment(s2,pi/2-tta,ltmp/cos(tta));
s4 = g.extendSegmentBySegment(s3,-pi/2+tta,dss-hs0-hs1);
s5 = g.extendSegmentBySegment(s4,-pi/2,wss/2);
s6 = g.extendSegmentBySegment(s5,pi/2,dsy);
s7 = g.extendSegmentBySegment(s6,pi/2,wss/2+wst/2);
s8 = g.extendSegmentBySegment(s7,pi/2,dss+dsy);

s9 = g.extendSegmentBySegment(s1,0,bs0/2);
s10 = g.extendSegmentBySegment(s9,pi/2,dss);

l1 = g.addLoop(s1,s2,s3,s4,s5,s6,s7,s8);
l2 = g.addRectangleLoop(wst/2+paperT,yShift+hs0+hs1+paperT,wss/2-paperT-bw0/2,dss-hs0-hs1-2*paperT);
l3 = g.addLoop(s9,s10,-s5,-s4,-s3,-s2);
g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3, l3, l2);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,255,137,39);
g.setFaceColor(name3,0,255,255);

% visualizations for debug
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
% stator tooth & coil

function emdlab_g2d_lib_tc9(g, ID, OD, Ns, Nc, wc, dc, bs0, hs0, tta, h0x, h0y, name1, name2, name3)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    ID (1,1) double {mustBePositive}
    OD (1,1) double {mustBePositive}
    Ns (1,1) double {mustBePositive,mustBeInteger}
    Nc (1,1) double {mustBePositive}
    wc (1,1) double {mustBePositive}
    dc (1,1) double {mustBePositive}
    bs0 (1,1) double {mustBePositive}
    hs0 (1,1) double {mustBePositive}
    tta (1,1) double {mustBePositive}
    h0x (1,1) double {mustBePositive}
    h0y (1,1) double {mustBePositive}
    name1 (1,:) char;
    name2 (1,:) char;
    name3 (1,:) char;

end

gamma_so = 2*asin(bs0/ID);
alpha_s = 2*pi/Ns;
tta = tta * pi/180;
ws = wc + 2*h0y;
ds = Nc*dc + (Nc+1)*h0x;
ltmp = ws/2-bs0/2;

x1 = (ID/2) * cos(gamma_so/2);
y1 = (ID/2) * sin(gamma_so/2);


e1 = g.addSegmentByCoordinates(x1,y1,x1+hs0,y1);
[e2,e2h] = g.extendSegmentBySegment(e1,pi/2-tta,ltmp/cos(tta));
e3 = g.extendSegmentBySegment(e2,-pi/2+tta,ds);
[e4, e4h] = g.extendSegmentBySegment(e3,-pi/2,ws/2);

e5 = g.addSegmentByCoordinates(e4h.ptr.p1.x,e4h.ptr.p1.y,OD/2,0);
e6 = g.extendSegmentByArc(e5,0,0,alpha_s/2);
e7 = g.extendArcBySegment(e6,pi/2,OD/2-ID/2);
e8 = g.extendSegmentByArc(e7,0,0,-alpha_s/2+gamma_so/2);

indexList = [];
xtmp = e2h.ptr.p1.x+ h0x;
for i = 1:Nc
    s1Index = g.addSegmentByCoordinates(xtmp,0,xtmp+dc,0);
    s2Index = g.addSegmentByCoordinates(xtmp+dc,0,xtmp+dc,wc/2);
    s3Index = g.addSegmentByCoordinates(xtmp+dc,wc/2,xtmp,wc/2);
    s4Index = g.addSegmentByCoordinates(xtmp,wc/2,xtmp,0);
    ltmp = g.addLoop(s1Index,s2Index,s3Index,s4Index);
    
    g.addFace(name2+string(i),ltmp);
    g.setFaceColor(name2+string(i),255,137,39);
    s5Index = g.addSegmentByCoordinates(xtmp+dc,0,xtmp+dc+h0x,0);
    indexList = [indexList,-s4Index,-s3Index,-s2Index,s5Index];
    xtmp = xtmp + dc + h0x;
end

l1 = g.addLoop(e1,e2,e3,e4,e5,e6,e7,e8);

e1tmp = g.addArcByCoordinates(0,0,x1,y1,ID/2,0,0);
e2tmp = g.addSegmentByCoordinates(ID/2,0,e2h.ptr.p1.x+ h0x,0);
indexList = [indexList,-e4,-e3,-e2,-e1,e1tmp,e2tmp];



l2 = g.addLoop(indexList);
% l3 = g.addLoop(e11,e1,e2,-e12,-e9);
% 
g.addFace(name1, l1);
% g.addFace(name2, l2);
g.addFace(name3, l2);
% 
g.setFaceColor(name1,200,200,200)

g.setFaceColor(name3,0,255,255)

end
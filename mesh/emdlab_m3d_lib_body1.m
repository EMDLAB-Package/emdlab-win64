function emdlab_m3d_lib_body1(m, pts1, pts2, p)

pts1 = emdlab_g2d_sortPointsCCW(pts1);

h0 = norm(pts1(2,:) - pts1(1,:));

g = emdlab_g2d_db;

ne = size(pts1,1) + 1;
pts1 = [0,0;pts1];
for i = 1:ne
    g.addPoint(pts1(i,1),pts1(i,2));
end
for i = 1:ne-1
    g.addSegment(i,i+1);
    g.edges(i).ptr.setNnodes(2);
end
g.addSegment(i+1,1);
g.edges(1).ptr.setMaxLength(h0);
g.edges(end).ptr.setMaxLength(h0);
g.addFace('z1', g.addLoop(1:ne));
g.setFaceColor('z1', 255, 0, 255);

tm = g.generateMesh('mg0');

w1 = p.Lstk/2;
w2 = w1 + p.housingEXT-p.hb;
dw = w2-w1;
z1 = linspace(w1, w2, ceil(dw/p.meshSizeZ));

w3 = w2 + p.hb;
dw = w3-w2;
z2 = linspace(w2, w3, ceil(dw/p.meshSizeZ));

w4 = w3 + +p.capT + p.shaftEXT;
dw = w4-w3;
z3 = linspace(w3, w4, ceil(dw/p.meshSizeZ));

tmp_a = atan(abs((p.Dsh-p.Db1)/2/(w2-w1)));
tmp_a = rad2deg(tmp_a);
if p.Dsh<p.Db1
    tmp_a = -tmp_a;
end



m.addMeshZone('shaft1', tm.mzs.z1.buildPrismMeshByExtrusion(p.z1));
[mzptr,llp] = tm.mzs.z1.buildPrismMeshByExtrusionZAxisDraft(z1,0,tmp_a);
m.addMeshZone('shaft2', mzptr);

tm.mzs.z1.nodes = llp(:,1:2);

m.addMeshZone('shaft3', tm.mzs.z1.buildPrismMeshByExtrusion(z2));
m.addMeshZone('shaft4', tm.mzs.z1.buildPrismMeshByExtrusion(z3));


pts2 = emdlab_g2d_sortPointsCCW(pts2);

h0 = norm(pts2(2,:) - pts2(1,:));

pts2 = flipud(pts2);

g = emdlab_g2d_db;

np = size(pts2,1);

u1 = pts2(1,:);
u1 = u1/norm(u1);
p1 = pts2(1,:) + p.housingT * u1;
s(1) = g.addSegmentByCoordinates(p1(1),p1(2),pts2(1,1),pts2(1,2));
g.edges(s(1)).ptr.setMaxLength(h0);
for i = 1:np-1
    s(end+1) = g.addSegmentByCoordinates(pts2(i,1),pts2(i,2),pts2(i+1,1),pts2(i+1,2));
    g.edges(s(end)).ptr.setNnodes(2);
end

u2 = pts2(end,:);
u2 = u2/norm(u2);
p2 = pts2(end,:) + p.housingT * u2;
s(end+1) = g.addSegmentByCoordinates(pts2(end,1),pts2(end,2),p2(1),p2(2));
g.edges(s(end)).ptr.setMaxLength(h0);

a1 = g.addArcByCoordinates(0,0,p2(1),p2(2),p1(1),p1(2));
g.edges(a1).ptr.setMaxLength(h0);

g.addFace('z1', g.addLoop([s,a1]));
g.setFaceColor('z1', 61, 145, 207);



tb = (p.Db2 - p.Db1)/2;
tc = p.capT;

p5 = g.addPoint(u2*(p.Db2/2+tc));
p6 = g.addPoint(u1*(p.Db2/2+tc));

p3 = g.addPoint(u2*p.Db2/2);
p4 = g.addPoint(u1*p.Db2/2);

p1 = g.addPoint(u2*(p.Db2/2-2*tb/3));
p2 = g.addPoint(u1*(p.Db2/2-2*tb/3));

p7 = g.addPoint(pts2(end,:));
p8 = g.addPoint(pts2(1,:));

o = g.addPoint(0,0);

a1 = g.addArc(o,p1,p2);
a2 = g.addArc(o,p3,p4);
a3 = g.addArc(o,p5,p6);

s1 = g.addSegment(p1,p3);
s2 = g.addSegment(p3,p5);
s3 = g.addSegment(p5,p7);

s4 = g.addSegment(p4,p2);
s5 = g.addSegment(p6,p4);
s6 = g.addSegment(p8,p6);

g.setMeshMaxLength(h0);

for i = 2:length(s)-1
    g.edges(s(i)).ptr.setNnodes(2);
end

g.addFace('z2', g.addLoop(s3,-fliplr(s(2:end-1)),s6,-a3));
g.addFace('z3', g.addLoop(s2,a3,s5,-a2));
g.addFace('z4', g.addLoop(s1,a2,s4,-a1));

tm = g.generateMesh('mg0');

z = linspace(p.Lstk/2, p.Lstk/2 + p.housingEXT, ceil(p.housingEXT/p.meshSizeZ));z(1)=[];
z1 = [p.z1,z];

m.addMeshZone('housing', tm.mzs.z1.buildPrismMeshByExtrusion1(z1));

tm.setMeshZoneColor('z3', 61, 145, 207);
tm.setMeshZoneColor('z4', 34, 177, 76);


z = linspace(p.Lstk/2 + p.housingEXT-p.hb, p.Lstk/2 + p.housingEXT, max(ceil(12/p.meshSizeZ),2));
m.addMeshZone('cap1', tm.mzs.z3.buildPrismMeshByExtrusion(z));

z = linspace(p.Lstk/2 + p.housingEXT-p.hb, p.Lstk/2 + p.housingEXT, max(ceil(12/p.meshSizeZ),2));
m.addMeshZone('b1', tm.mzs.z4.buildPrismMeshByExtrusion(z));

tm.aux_unify('z');

z = linspace(z1(end), p.Lstk/2 + p.housingEXT+tc, max(ceil(tc/p.meshSizeZ),2));

m.addMeshZone('cap2', tm.mzs.z.buildPrismMeshByExtrusion1(z));

pts3 = g.edges(a1).ptr.getMeshNodes;


g = emdlab_g2d_db;
pts1 = pts1(2:end,:);

u = pts1;
u = u./vecnorm(u,2,2);
pts1 = pts1 - u * (p.Dsh/2 - p.Db1/2);
n1 = size(pts1,1);
n2 = size(pts3,1);

pts = [pts3;flipud(pts1)];
g.addClosedPolylineLoop(pts(:,1),pts(:,2));
g.addFace('z1', 1);

for i = 1:(n1+n2)
    g.edges(i).ptr.setNnodes(2);
end
g.edges(n2).ptr.setMaxLength(h0);
g.edges(n1+n2).ptr.setMaxLength(h0);

tm = g.generateMesh('mg0');

tm.setMeshZoneColor('z1', 34, 177, 76);

z = linspace(w2, w3, max(ceil((w3-w2)/p.meshSizeZ),2));

m.addMeshZone('b2', tm.mzs.z1.buildPrismMeshByExtrusion(z));

m.joinMeshZones('bearing','b1','b2')
m.joinMeshZones('cap', 'cap1', 'cap2')
m.joinMeshZones('shaft', 'shaft1', 'shaft2', 'shaft3', 'shaft4')

end
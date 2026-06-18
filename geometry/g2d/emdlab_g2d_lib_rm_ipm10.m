% rotor & magnet
% outer rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_rm_ipm1(g, ID, OD, poles, Kair, g_wry, g_dm, g_wm, wrribs, wtrib, name1, name2, name3)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    ID (1,1) double {mustBePositive}
    OD (1,1) double {mustBePositive}
    poles (1,1) double {mustBePositive,mustBeInteger}
    Kair (1,1) double {mustBePositive}
    g_wry (1,:) double {mustBePositive}
    g_dm (1,:) double {mustBePositive}
    g_wm (1,:) double {mustBePositive}
    wrribs (1,:) double {mustBePositive}
    wtrib (1,1) double {mustBePositive}
    name1 (1,:) char;
    name2 (1,:) char;
    name3 (1,:) char;

end

% the number of flux barrier layers
Nfb = length(g_wry);
% if length(g_dm) ~= (Nfb-1)
%     error('length of g_dm must be Nfb-1.');
% end

alpha_p = 2*pi/poles;

lrIndex = g.addAnnularSectorLoop(ID/2, OD/2, 0, alpha_p/2, 0, 0);

o = emdlab_g2d_point(0,0);

% **********************************************************************
% calculate contour lines
tmp = 1;
for i = 1:numel(g_wry)
    tmp = tmp + prod(g_wry(1:i));
end
gv_wry = zeros(1,Nfb+1);
gv_wry(1) = (OD/2-ID/2)*(1-Kair)/tmp;
for i = 2:Nfb+1
    gv_wry(i) = gv_wry(i-1)*g_wry(i-1);
end

% total thickness of bagv_riers
dmt = (OD/2-ID/2)*Kair;
tmp = 1;
for i = 1:numel(g_dm)
    tmp = tmp + prod(g_dm(1:i));
end
gv_dm = zeros(1,Nfb);
gv_dm(1) = dmt / tmp;
for i = 2:Nfb
    gv_dm(i) = gv_dm(i-1)*g_dm(i-1);
end

% maximum magnet width
wm_max = tan(alpha_p/2) * (ID/2+gv_wry(1)+gv_dm(1))*2;

wm = zeros(1,Nfb);
wm(1) = g_wm(1)*wm_max;
for i = 2:Nfb
    wm(i) = wm(i-1)*g_wm(i);
end

ux = cos(alpha_p/2);
uy = sin(alpha_p/2);

for i = 1:Nfb

    x1 = (ID/2) + sum(gv_wry(1:i)) + sum(gv_dm(1:i-1));
    x2 = x1 + gv_dm(i);

    y1 = wrribs(i)/2;
    y2 = wrribs(i)/2;

    x3 = x2;
    y3 = y2+wm(i)/2;

    

    s1Index = g.addSegmentByCoordinates(x1,y1,x2,y2);
    s2Index = g.addSegmentByCoordinates(x2,y2,x3,y3);

    s2Index = g.splitEdge(s2Index,[0.1,0.8]);

    sm1Index = g.extendSegmentBySegment(s2Index(1),pi/2,gv_dm(i));
    sm2Index = g.extendSegmentBySegment(s2Index(2),pi/2,gv_dm(i));

    err_f = @(t) norm([x3,y3]+[ux,uy]*t)-OD/2+wtrib;
    t = fzero(err_f,OD/2);

    x4 = x3 +t*ux;
    y4 = y3 +t*uy;

    s3Index = g.addSegmentByCoordinates(x3,y3,x4,y4);

    s3Index = g.splitEdge(s3Index,[0.2,0.6]);

    sm3Index = g.extendSegmentBySegment(s3Index(1),pi/2,gv_dm(i));
    sm4Index = g.extendSegmentBySegment(s3Index(2),pi/2,gv_dm(i));
    



    x6 = x1;
    y6 = y1+wm(i)/2+gv_dm(i)/tan((pi/2+alpha_p/2)/2);

    err_f = @(t) norm([x6,y6]+[ux,uy]*t)-OD/2+wtrib;
    t = fzero(err_f,OD/2);

    x5 = x6 +t*ux;
    y5 = y6 +t*uy;

    s4Index = g.addArcByCoordinates(0,0,x4,y4,x5,y5,1);

    s5Index(1) = g.addSegmentByCoordinates(x5,y5,g.edges(sm4Index).ptr.p1.x,g.edges(sm4Index).ptr.p1.y);
    s5Index(2) = g.addSegmentByCoordinates(g.edges(sm4Index).ptr.p1.x,g.edges(sm4Index).ptr.p1.y,...
        g.edges(sm3Index).ptr.p1.x,g.edges(sm3Index).ptr.p1.y);
    s5Index(3) = g.addSegmentByCoordinates(g.edges(sm3Index).ptr.p1.x,g.edges(sm3Index).ptr.p1.y,x6,y6);


    s6Index(1) = g.addSegmentByCoordinates(x6,y6,g.edges(sm2Index).ptr.p1.x,g.edges(sm2Index).ptr.p1.y);
    s6Index(2) = g.addSegmentByCoordinates(g.edges(sm2Index).ptr.p1.x,g.edges(sm2Index).ptr.p1.y,...
        g.edges(sm1Index).ptr.p1.x,g.edges(sm1Index).ptr.p1.y);
    s6Index(3) = g.addSegmentByCoordinates(g.edges(sm1Index).ptr.p1.x,g.edges(sm1Index).ptr.p1.y,x1,y1);


   lair1 = g.addLoop(s1Index,s2Index(1),sm1Index, s6Index(3));
   lair2 = g.addLoop(s2Index(3),s3Index(1),sm3Index, s5Index(3), s6Index(1),-sm2Index);
   lair3 = g.addLoop(s3Index(3),s4Index,s5Index(1),-sm4Index);

   g.addFace(name3+string(3*i-2), lair1);
   g.addFace(name3+string(3*i-1), lair2);
   g.addFace(name3+string(3*i), lair3);

   g.setFaceColor(name3+string(3*i-2),0,255,255)
   g.setFaceColor(name3+string(3*i-1),0,255,255)
   g.setFaceColor(name3+string(3*i),0,255,255)

   lm1 = g.addLoop(s2Index(2),sm2Index, s6Index(2),-sm1Index);
   lm2 = g.addLoop(s3Index(2),sm4Index, s5Index(2),-sm3Index);

   g.addFace(name2+"c"+string(i), lm1);
   g.addFace(name2+"s"+string(i), lm2);
   g.setFaceColor(name2+"c"+string(i),170,218,24)
    g.setFaceColor(name2+"s"+string(i),170,218,24)

   l111 = g.addLoop(s1Index,s2Index,s3Index,s4Index,s5Index,s6Index);
   lrIndex(end+1) = l111;

end


g.addFace(name1, lrIndex);
% g.addFace(name2, l2);
% g.addFace(name3, l3);
% 
% set mesh zone colors
g.setFaceColor(name1,200,200,200)
% g.setFaceColor(name2,28,255,28)
% g.setFaceColor(name3,0,255,255)

end
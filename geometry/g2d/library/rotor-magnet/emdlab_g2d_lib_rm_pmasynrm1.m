% rotor & magnets & flux barriers
% permanent magnet assissted synchronous reluctance motor

function emdlab_g2d_lib_rm_pmasynrm1(g, ID, OD, p, Kair, g_wry, g_dm, g_wm, wrribs, wtrib, name1, name2, name3)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    ID (1,1) double {mustBePositive}
    OD (1,1) double {mustBePositive}
    p (1,1) double {mustBePositive,mustBeInteger}
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
gv_Nfb = length(g_wry);
if length(g_dm) ~= (gv_Nfb-1)
    error('length of g_dm must be ')
end

% parameters
gv_Bav = 0.45;
gv_ri = 24/2;
gv_ro = 74/2;
gv_p = 4;
gv_Kair = 0.4;
gv_wtrib = 0.6;
gv_Nfb = 4;
gv_gwry = 0.9*ones(1,gv_Nfb);
gv_gdb = 0.7*ones(1,gv_Nfb - 1);

% creation of the geometry
g = emdlab_g2d_db;
g.addAnnularSectorLoop(gv_ri, gv_ro, 0, pi/gv_p, 0, 0);
g.addFace('rotor', 1);
g.setMeshMaxLength(1);

% mesh generation
m = g.generateMesh('gmsh');
% m.addMaterial('es_M400_50A');
% m.setMaterial('rotor', 'es_M400_50A');

% define solver
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');

% set boundary conditions
k1 = m.getfbnioc([0,0], gv_ri);
k1 = [k1; m.getfbniol_p0p1([0,0], [cos(pi/gv_p),sin(pi/gv_p)])];
k2 = m.getfbnioc([0,0], gv_ro);

% inner surface
s.setAzBC(k1,0);

% outer surface
p = m.nodes(k2,:);
p = atan2(p(:,2),p(:,1));
s.setAzBC(k2,(-2/gv_p)*(gv_ro/1000)*(pi*gv_Bav/2)*cos((gv_p/2)*p));

% run solver
s.setMonitor(0);
s.solve;
s.plotBmagF(2*gv_Nfb+2);
clim([0,pi*gv_Bav/2]);

% **********************************************************************
% calculate contour lines
tmp = 1;
for i = 1:numel(gv_gwry)
    tmp = tmp + prod(gv_gwry(1:i));
end
gv_wry = zeros(1,gv_Nfb+1);
gv_wry(1) = (gv_ro-gv_ri)*(1-gv_Kair)/tmp;
for i = 2:gv_Nfb+1
    gv_wry(i) = gv_wry(i-1)*gv_gwry(i-1);
end

% total thickness of bagv_riers
dbt = (gv_ro-gv_ri)*gv_Kair;
tmp = 1;
for i = 1:numel(gv_gdb)
    tmp = tmp + prod(gv_gdb(1:i));
end
gv_db = zeros(1,gv_Nfb);
gv_db(1) = dbt / tmp;
for i = 2:gv_Nfb
    gv_db(i) = gv_db(i-1)*gv_gdb(i-1);
end

% calculation of starting points
xStart = zeros(1,2*gv_Nfb);
yStart = zeros(1,2*gv_Nfb);

for i = 1:gv_Nfb

    xStart(2*i-1) = gv_ri + sum(gv_wry(1:i)) + sum(gv_db(1:i-1));
    xStart(2*i) = xStart(2*i-1) + gv_db(i);

end

% define cell variables to store contour line coordinates
xb = cell(1,2*gv_Nfb);
yb = cell(1,2*gv_Nfb);

alpha_tmp = linspace(0,pi,1000);

% calculation loop for all contour lines
for i = 1:2*gv_Nfb

    xb{i}(1) = xStart(i);
    yb{i}(1) = yStart(i);

    % calculation of Az at initial point
    Ac = s.getAzOnPoints(xStart(i), yStart(i));

    % calculation loop for one contour line
    while true

        xp = (gv_wtrib*0.9)*cos(alpha_tmp) + xb{i}(end);
        yp = (gv_wtrib*0.9)*sin(alpha_tmp) + yb{i}(end);

        Ap = s.getAzOnPoints(xp, yp);

        err = abs(Ap-Ac);
        [~, index] = min(err);

        xb{i}(end+1) = xp(index);
        yb{i}(end+1) = yp(index);

        % condition to find end point
        if norm([xp(index),yp(index)]) > gv_ro-gv_wtrib

            t_tmp = linspace(0,pi/gv_p,100000);
            xp = (gv_ro-gv_wtrib)*cos(t_tmp);
            yp = (gv_ro-gv_wtrib)*sin(t_tmp);

            Ap = s.getAzOnPoints(xp, yp);

            err = abs(Ap-Ac);
            [~, index] = min(err);

            xb{i}(end) = xp(index);
            yb{i}(end) = yp(index);
            break;

        end

    end

end

figure;
axis off equal;
hold on;
plot(xStart,yStart,'*')
for i = 1:2*gv_Nfb
    plot(xb{i},yb{i},'k', 'LineWidth',2);    
end

g2 = emdlab_g2d_db;
for i = 1:gv_Nfb
    g2.addSegmentByCoordinates(xb{2*i-1}(1),0,xb{2*i}(1),0);
    g2.addSplineByCoordinates(xb{2*i},yb{2*i});
    g2.addArcByCoordinates(0,0,xb{2*i}(end),yb{2*i}(end),xb{2*i-1}(end),yb{2*i-1}(end),1);
    g2.addSplineByCoordinates(fliplr(xb{2*i-1}),fliplr(yb{2*i-1}));
    g2.addLoop(4*(i-1) + (1:4));
    g2.addFace(['fb',num2str(i)], i);
end

g2.addSegmentByCoordinates(gv_ri,0,xb{1}(1),0);
for i = 1:gv_Nfb-1
    g2.addSegmentByCoordinates(xb{2*i}(1),0,xb{2*i+1}(1),0);
end
e = g2.addSegmentByCoordinates(xb{end}(1),0,gv_ro,0);
e = g2.extendSegmentByArc(e,0,0,pi/gv_p);
e = g2.extendArcBySegment(e,pi/2,gv_ro-gv_ri);
g2.extendSegmentByArc(e,0,0,-pi/gv_p);

tmp = 1:24;
tmp = setdiff(tmp,1:4:13);
tmp(1:12) = -tmp(1:12);
g2.addLoop(17,-4,-3,-2,18,-8,-7,-6,19,-12,-11,-10,20,-16,-15,-14,21:24);
g2.addFace('rotor',5);
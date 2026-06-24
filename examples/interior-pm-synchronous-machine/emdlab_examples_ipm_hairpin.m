% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design parameters
gv_ISD = 95;
gv_OSD = 150;
gv_g = 0.7;
gv_Lstk = 100;
gv_Dsh = 40;
gv_Ns = 54;
gv_Nc = 6;
gv_p = 6;
gv_Hc = -922100;

% dependent variables
alpha_p = 2*pi/gv_p;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc9(g,gv_ISD,gv_OSD,gv_Ns,gv_Nc,2.8,1.8,1.5,0.8,15,0.3,0.3,'stator','sc','sap');
emdlab_g2d_lib_rm_ipm13(g, gv_Dsh, gv_ISD-2*gv_g, gv_p, 0.25, [0.8,1], [1]*.8, [0.6,0.7], [1,1]*0.5, 0.6, 'rotor', 'magnet', 'rap');

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [1.5,gv_g,2], r, 'linear', 'extrap');
g.setMeshLengthByRadialFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');

% add materials
m.addMaterial('m330', emdlab_mlib_es_M300_35A);
m.addMaterial('copper', emdlab_mlib_copper);

% set materials
m.setMaterial('rotor','m330');
m.setMaterial('stator','m330');

% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns);
m.aux_cmxjcrj('sap',gv_Ns);
for i = 1:gv_Nc
m.aux_cmxjcrj('sc'+string(i),gv_Ns);
end
m.aux_cmxjcrj('rotor',gv_p);
m.aux_unify('rap');
m.aux_cmxjcrj('rap',gv_p);
m.joinMeshZones('magnetc', 'magnetc'+string(1:2))
m.aux_cmxjcr('magnetc',gv_p);
m.joinMeshZones('magnets', 'magnets'+string(1:2))
m.aux_cmxcr('magnets',gv_p);

% add air gap mesh
m.aux_addCircularAirGap('ag', 0,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,2);

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);
magVector = [1,0];
u1 = [cos(alpha_p/2),sin(alpha_p/2)];
u1s = emdlab_g2d_rotatePoints(u1,-pi/2);
u2s = emdlab_g2d_rotatePoints([u1(1),-u1(2)],pi/2);

% set magnetization
for i = 1:2:gv_p
    magVector1 = emdlab_g2d_rotatePoints(magVector,(i-1)*alpha_p);
    magVector2 = emdlab_g2d_rotatePoints(u1s,(i-1)*alpha_p);
    magVector3 = emdlab_g2d_rotatePoints(u2s,(i-1)*alpha_p);

    s.setMagnetization(['magnetc',num2str(i)],gv_Hc,magVector1);
    s.setMagnetization(['magnets1',num2str(i)],gv_Hc,magVector2);
    s.setMagnetization(['magnets2',num2str(i)],gv_Hc,magVector3);

    m.setMeshZoneColor(['magnetc',num2str(i)],0,125,223);
    m.setMeshZoneColor(['magnets1',num2str(i)],0,125,223);
    m.setMeshZoneColor(['magnets2',num2str(i)],0,125,223);
end
for i = 2:2:gv_p
    magVector1 = emdlab_g2d_rotatePoints(magVector,(i-1)*alpha_p);
    magVector2 = emdlab_g2d_rotatePoints(u1s,(i-1)*alpha_p);
    magVector3 = emdlab_g2d_rotatePoints(u2s,(i-1)*alpha_p);

    s.setMagnetization(['magnetc',num2str(i)],gv_Hc,-magVector1);
    s.setMagnetization(['magnets1',num2str(i)],gv_Hc,-magVector2);
    s.setMagnetization(['magnets2',num2str(i)],gv_Hc,-magVector3);

    m.setMeshZoneColor(['magnetc',num2str(i)],255,70,70);
    m.setMeshZoneColor(['magnets1',num2str(i)],255,70,70);
    m.setMeshZoneColor(['magnets2',num2str(i)],255,70,70);
end

% apply boundary conditions: zero vector potential on all boundary nodes
s.setAzBC(m.getfbn, 0);

% solve and plot results
s.solve(20);
g.showSketch;
m.showg;
m.showmzs;
f = s.plotBmagF(14, 'rotor', 'stator');
s.plotBrBtOnCircle(0, 0, gv_ISD/2-gv_g/2, 1000);

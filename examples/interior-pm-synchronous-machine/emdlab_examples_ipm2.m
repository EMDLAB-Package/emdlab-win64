%{
note: calculation of the rotor field of an v-shape ipm motor
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 84;
gv_OSD = 140;
gv_Lstk = 60;
gv_Ns = 12;
gv_p = 8;
gv_wst = 12;
gv_dss = 20;
gv_bs0 = 5;
gv_hs0 = 2;
gv_tta = 15;
gv_bw0 = 1;
gv_Dsh = 45;
gv_dm = 4;
gv_gap = 0.5;
gv_wtrib = 1;
gv_wrrib = 0;
gv_bm0 = 5;
gv_alpha_v = 140;
gv_Hc = -922100;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc21(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, gv_bw0, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_ipm5(g,gv_Dsh,gv_ISD-2*gv_gap,gv_p,gv_dm,gv_alpha_v,gv_wtrib,gv_wrrib,gv_bm0,'rotor','magnet','rap')

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [3,gv_gap/2,3], r, 'linear','extrap');
g.setMeshLengthByRadialFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');

% add materials
m.addMaterial('m330', emdlab_mlib_es_M330_35A);
m.addMaterial('copper', emdlab_mlib_copper);

% set materials
m.setMaterial('rotor','m330');
m.setMaterial('stator','m330');
m.setMaterial('sc','copper');

% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns)
m.aux_cmxjcrj('sap',gv_Ns)
m.aux_cmxjcrj('rotor',gv_p)
m.aux_unify('rap');
m.aux_cmxjcrj('rap',gv_p)
m.aux_cmxjcrj('sc',gv_Ns)
m.aux_cmxcr('magnet',gv_p)

% generate air gap mesh
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_gap,0,0,gv_ISD/2,2)

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set magnetization
tmp = pi/2-deg2rad(gv_alpha_v/2);
u1 = [cos(tmp),-sin(tmp)];
u2 = [u1(1),-u1(2)];
for i = 1:2:gv_p
    u11 = emdlab_g2d_rotatePoints(u1,(i-1)*2*pi/gv_p);
    s.setMagnetization(['magnet1',num2str(i)],gv_Hc,u11);
    m.setMeshZoneColor(['magnet1',num2str(i)],0,125,223);
    u22= emdlab_g2d_rotatePoints(u2,(i-1)*2*pi/gv_p);
    s.setMagnetization(['magnet2',num2str(i)],gv_Hc,u22);
    m.setMeshZoneColor(['magnet2',num2str(i)],0,125,223);
end
u1 = -emdlab_g2d_rotatePoints(u1,2*pi/gv_p);
u2 = -emdlab_g2d_rotatePoints(u2,2*pi/gv_p);
for i = 2:2:gv_p
    u11 = emdlab_g2d_rotatePoints(u1,(i-2)*2*pi/gv_p);
    s.setMagnetization(['magnet1',num2str(i)],gv_Hc,u11);
    m.setMeshZoneColor(['magnet1',num2str(i)],255,70,70);
    u22= emdlab_g2d_rotatePoints(u2,(i-2)*2*pi/gv_p);
    s.setMagnetization(['magnet2',num2str(i)],gv_Hc,u22);
    m.setMeshZoneColor(['magnet2',num2str(i)],255,70,70);
end

% apply boundary conditions: zero vector potential on all boundary nodes
s.setAzBC(m.getfbn, 0);

% run solver
s.solve(20);

% visualize solution
s.plotBmagF(16,'rotor','stator');
m.showg;
s.plotBrBtOnCircle(0, 0, gv_ISD/2-gv_gap/2, 1000);
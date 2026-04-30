%{
note: calculation of the rotor field of an v-shape ipm motor
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 45;
gv_OSD = 75;
gv_Lstk = 60;
gv_Ns = 24;
gv_p = 4;
gv_wst = 2.6;
gv_dss = 10;
gv_bs0 = 1.5;
gv_hs0 = 0.6;
gv_tta = 25;
gv_Dsh = 14;
gv_dm = 3;
gv_g = 0.5;
gv_embrace = 0.85;
gv_ra = 20;
gv_Hc = -922100;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc1(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_ipm8(g,gv_Dsh,gv_ISD-2*gv_g,gv_p,gv_embrace,gv_ra,gv_dm,0.8,0.6,0.6,'rotor','magnet','rap')

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [2,gv_g/2,2], r, 'linear','extrap');
g.setMeshLengthByRadialFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');
m.joinMeshZones('rap','rap1','rap2')

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
m.aux_cmxjcrj('rap',gv_p)
m.aux_cmxjcr('sc',gv_Ns)
m.aux_cmxcr('magnet',gv_p)

% generate air gap mesh
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,2)

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set magnetization
tmp = gv_ra*pi/180;
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
s.plotBrBtOnCircle(0, 0, gv_ISD/2-gv_g/2, 1000);
s.plotBmagF;
clim([0,1.8])
m.showm
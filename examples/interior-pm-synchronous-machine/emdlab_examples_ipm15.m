%{
note: calculation of the rotor field of an u-shape ipm motor
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 85;
gv_OSD = 140;
gv_Lstk = 60;
gv_gap = 0.7;
gv_Ns = 36;
gv_p = 6;
gv_wst = 4;
gv_dss = 15;
gv_bs0 = 1.5;
gv_hs0 = 1;
gv_tta = 15;
gv_bw0 = 1;
gv_Dsh = 35;
gv_dm = 3.5;
gv_wsm = 8;
gv_w1 = 8;
gv_w2 = 0.7;
gv_w3 = 1;
gv_w4 = 1;
gv_Hc = -922100;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_ipm15(g,gv_Dsh,gv_ISD-2*gv_gap,gv_p,gv_dm,gv_wsm,gv_w1,gv_w2,gv_w3,gv_w4,'rotor','magnet','rap')

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [2,gv_gap/2,2], r, 'linear','extrap');
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
m.joinMeshZones('rap','rap1','rap2')
m.aux_cmxjcrj('rap',gv_p)
m.aux_cmxjcrj('sc',gv_Ns)
m.changeMeshZoneName('magnet1', 'magnet');
m.changeMeshZoneName('magnet2', 'smagnet');
m.aux_cmxjcr('magnet',gv_p)
m.aux_cmxcr('smagnet',gv_p)

% generate air gap mesh
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_gap,0,0,gv_ISD/2,2)

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set magnetization
tmp = pi/2-pi/gv_p;
u1 = [cos(tmp),-sin(tmp)];
u2 = [u1(1),-u1(2)];
for i = 1:2:gv_p
    u11 = emdlab_g2d_rotatePoints(u1,(i-1)*2*pi/gv_p);
    s.setMagnetization(['smagnet1',num2str(i)],gv_Hc,u11);
    m.setMeshZoneColor(['smagnet1',num2str(i)],0,125,223);
    u22= emdlab_g2d_rotatePoints(u2,(i-1)*2*pi/gv_p);
    s.setMagnetization(['smagnet2',num2str(i)],gv_Hc,u22);
    m.setMeshZoneColor(['smagnet2',num2str(i)],0,125,223);
    uMag = emdlab_g2d_rotatePoints([1,0],(i-1)*2*pi/gv_p);
    s.setMagnetization(['magnet',num2str(i)],gv_Hc,uMag);
    m.setMeshZoneColor(['magnet',num2str(i)],0,125,223);
end
u1 = -emdlab_g2d_rotatePoints(u1,2*pi/gv_p);
u2 = -emdlab_g2d_rotatePoints(u2,2*pi/gv_p);
for i = 2:2:gv_p
    u11 = emdlab_g2d_rotatePoints(u1,(i-2)*2*pi/gv_p);
    s.setMagnetization(['smagnet1',num2str(i)],gv_Hc,u11);
    m.setMeshZoneColor(['smagnet1',num2str(i)],255,70,70);
    u22= emdlab_g2d_rotatePoints(u2,(i-2)*2*pi/gv_p);
    s.setMagnetization(['smagnet2',num2str(i)],gv_Hc,u22);
    m.setMeshZoneColor(['smagnet2',num2str(i)],255,70,70);
    uMag = emdlab_g2d_rotatePoints([1,0],(i-1)*2*pi/gv_p);
    s.setMagnetization(['magnet',num2str(i)],gv_Hc,-uMag);
    m.setMeshZoneColor(['magnet',num2str(i)],255,70,70);
end

% apply boundary conditions: zero vector potential on all boundary nodes
s.setAzBC(m.getfbn, 0);

% run solver
s.solve(20);

% visualize solution
s.plotBmagF(16,'rotor','stator');
m.showg;
s.plotBrBtOnCircle(0, 0, gv_ISD/2-gv_gap/2, 1000);
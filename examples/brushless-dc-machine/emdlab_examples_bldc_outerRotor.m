%{
note: rotor field calculation for an outer-rotor BLDC motor
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% dimensions & parameters
gv_ISD = 50; % inner stator diamater [mm]
gv_OSD = 85; % outer stator diameter [mm]
gv_Lstk = 25; % stack length [mm]
gv_Ns = 24; % number of stator slots
gv_p = 20; % number of rotor poles
gv_wst = 4.5; % width of stator tooth [mm]
gv_dss = 10; % depth of stator slot [mm]
gv_bs0 = 3; % stator slot opening [mm] 
gv_hs0 = 1; % stator slot opening height [mm]
gv_tta = 30; % tooth tip angle of stator slot [deg]
gv_ORD = 96; % outer rotor diameter [mm]
gv_dm = 1.5; % magnet depth [mm]
gv_g = 0.7; % airgap thickness [mm]
gv_embrace = 0.89; % magnet embrace ratio -> lower than one
gv_Hc = -922100; % magnet coercive force [A/m]

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc2(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_spm2(g, gv_OSD+2*gv_g, gv_ORD, gv_p, gv_dm, gv_embrace, 'rotor', 'magnet', 'rap');

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_ISD/2,gv_OSD/2,gv_OSD/2+gv_g+gv_dm,gv_ORD/2], [1.5,0.5,0.5,1], r, 'linear','extrap');
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
m.aux_cmxjcrj('rap',gv_p)
m.aux_cmcr('sc',[cos(2*pi/gv_Ns),sin(2*pi/gv_Ns)],gv_Ns)
m.aux_cmxjcr('magnet',gv_p)

% add circular air gap
m.aux_addCircularAirGap('ag',0,0,gv_OSD/2,0,0,gv_OSD/2+gv_g,2);

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set magnetization
for i = 1:2:gv_p
s.setMagnetization(['magnet',num2str(i)],gv_Hc,'r');
m.setMeshZoneColor(['magnet',num2str(i)],0,125,223);
end
for i = 2:2:gv_p
s.setMagnetization(['magnet',num2str(i)],gv_Hc,'-r');
m.setMeshZoneColor(['magnet',num2str(i)],255,70,70);
end

% apply boundary conditions: zero vector potential on all boundary nodes
s.setAzBC(m.getfbn, 0);

% run solver
s.solve;

% visualize solution
s.gui;
s.plotBrBtOnCircle(0, 0, gv_OSD/2+gv_g/2, 1000);



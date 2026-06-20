%{
note: calculation of the rotor field of a single-rotor single-stator
axial-flux permanent magnet motor
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 150;
gv_OSD = 250;
gv_Ns = 18;
gv_p = 20;
gv_wss = 14;
gv_dss = 40;
gv_dsy = 10;
gv_bs0 = 5;
gv_hs0 = 1.5;
gv_tta = 25;
gv_Ntc = 85;
gv_g = 1;
gv_Iph = 1.9;
gv_Nfb = 4;
gv_paperT = 1;
gv_bw0 = 2;
gv_dm = 6;
gv_dry = 6;
gv_embrace = 0.8;
gv_Hc = -922100;

% dependent variables
gv_Lstk = (gv_OSD-gv_ISD)/2;
gv_D = (gv_OSD+gv_ISD)/2;
gv_L = pi*gv_D;
gv_taus = gv_L/gv_Ns;
gv_wst = gv_taus - gv_wss;
gv_taup = gv_L/gv_p;

% define geometry data base
g = emdlab_g2d_db;

% add geometries from templates
emdlab_g2d_lib_tc51(g,gv_wst,gv_wss,gv_dsy,gv_dss,gv_bs0,gv_hs0,gv_tta,0,gv_g/2+gv_dss+gv_dsy,1,'stator','sc','sap');
emdlab_g2d_lib_rm_spm51(g,-gv_g/2,gv_dm,gv_dry,gv_taup,gv_embrace,'rotor','magnet','rap');

% setting the wireframe mesh by mesh size function
f_mesh = @(y) interp1([-gv_g-gv_dm-gv_dry,0,gv_g+gv_dss+gv_dsy],[gv_dry/2,gv_g,gv_dsy/2],y, 'linear','extrap');
g.setMeshLengthByYFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');

% add materials
m.addMaterial('m530', emdlab_mlib_es_M530_50A);
m.addMaterial('copper', emdlab_mlib_copper);

% set mesh zone materials
m.setMaterial('rotor','m530');
m.setMaterial('stator','m530');
m.setMaterial('sc','copper');

% generate full mesh
m.aux_cmyjcshj('stator',gv_Ns,gv_taus,0);
m.aux_cmyjcshj('sap',gv_Ns,gv_taus,0);
m.aux_cmycsh('sc',gv_Ns,gv_taus,0);
m.aux_cmyjcshj('rotor',gv_p,gv_taup,0);
m.aux_cmyjcshj('rap',gv_p,gv_taup,0);
m.aux_cmyjcsh('magnet',gv_p,gv_taup,0);

% add airgap mesh
m.aux_addLineAirGap('ag',0,0,1,0,gv_g,2)

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set magnetization
for i = 1:2:gv_p
s.setMagnetization(['magnet',num2str(i)],gv_Hc,[0,1]);
m.setMeshZoneColor(['magnet',num2str(i)],0,125,223);
end
for i = 2:2:gv_p
s.setMagnetization(['magnet',num2str(i)],gv_Hc,[0,-1]);
m.setMeshZoneColor(['magnet',num2str(i)],255,70,70);
end

% apply boundary conditions
k1 = m.getfbniol_p0x([0,-gv_g/2-gv_dm-gv_dry]);
k2 = m.getfbniol_p0x([0,gv_g/2+gv_dss+gv_dsy]);
s.setAzBC([k1;k2],0);
[km,ks] = m.splitShift(setdiff(m.getfbn,[k1;k2]), [gv_L,0]);
s.setEvenPeriodicBC(km,ks);

% run solver
s.solve;

% visualize solution
s.gui;



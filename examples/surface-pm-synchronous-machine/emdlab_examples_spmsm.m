%{
note: calculation of the rotor field of a surface-mounted PMSM
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 74;
gv_OSD = 125;
gv_Lstk = 70;
gv_Ns = 36;
gv_p = 4;
gv_wst = 3;
gv_dss = 15;
gv_bs0 = 2.4;
gv_hs0 = 0.6;
gv_tta = 25;
gv_Dsh = 28;
gv_dm = 3;
gv_g = 1;
gv_embrace = 0.89;
gv_Hc = -922100;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc1(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_spm1(g, gv_Dsh, gv_ISD-2*gv_g, gv_p, gv_dm, gv_embrace, 'rotor', 'magnet', 'rap')

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_dm-gv_g,gv_ISD/2,gv_OSD/2], [4,0.5,0.5,4], r, 'linear','extrap');
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
m.aux_cmxjcr('sc',gv_Ns)
m.aux_cmxjcr('magnet',gv_p)

% generate air gap mesh
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,2)

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
m.showmzs;
f = s.plotBmagF(14, 'rotor', 'stator', 'magnet' + string(1:gv_p));
s.plotBrBtOnCircle(0, 0, gv_ISD/2-gv_g/2, 1001);

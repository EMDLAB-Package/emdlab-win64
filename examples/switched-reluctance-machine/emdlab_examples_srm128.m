%{
note: magneto-static analysis of a 12/8 SRM by exciting one stator phase 
for a specific rotor position
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% dimensions & parameters
gv_Ns = 12;
gv_Nr = 8;
gv_Dsh = 24;
gv_ISD = 70;
gv_OSD = 125;
gv_Lstk = 70;
gv_g = 0.45;
gv_gbetas = 0.5;
gv_gbetar = 0.38;
gv_gwsy = 1;
gv_gwry = 0.9;
gv_Iphase = 3;
gv_Ntc = 200;
rotorPosition = -10;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc_srm1(g, gv_ISD, gv_OSD, gv_Ns, gv_gbetas, gv_gwsy, 'stator', 'sca');
emdlab_g2d_lib_tc_srm2(g, gv_Dsh, gv_ISD-2*gv_g, gv_Nr, gv_gbetar, gv_gwry, 'rotor', 'rap');

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [5,2,5], r, 'linear','extrap');
g.setMeshLengthByRadialFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');
m.setMeshZoneColor('rap', 0, 255, 255)

% add materials
m.addMaterial('m530', emdlab_mlib_es_M530_50A);
m.addMaterial('copper', emdlab_mlib_copper);

% set materials
m.setMaterial('rotor','m530');
m.setMaterial('stator','m530');
m.setMaterial('sca','copper');

% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns)
m.aux_cmxjcrj('rotor',gv_Nr)
m.aux_cmxjcrj('rap', gv_Nr);
m.aux_cmxcr('sca', gv_Ns);

% rotate moving mesh zones
m.rotateMeshZone('rotor', rotorPosition*pi/180);
m.rotateMeshZone('rap', rotorPosition*pi/180);

% add airgap mesh
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,2);

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% define new winding
s.defineCoil('phaseA');
s.addMeshZone2Coil('phaseA','sca11',gv_Ntc,1);
s.addMeshZone2Coil('phaseA','sca21',gv_Ntc,-1);
s.addMeshZone2Coil('phaseA','sca14',gv_Ntc,-1);
s.addMeshZone2Coil('phaseA','sca24',gv_Ntc,1);
s.addMeshZone2Coil('phaseA','sca17',gv_Ntc,1);
s.addMeshZone2Coil('phaseA','sca27',gv_Ntc,-1);
s.addMeshZone2Coil('phaseA','sca110',gv_Ntc,-1);
s.addMeshZone2Coil('phaseA','sca210',gv_Ntc,1);

% set winding current
s.setCoilCurrent('phaseA', gv_Iphase);

% apply boundary conditions
s.setAzBC(m.getfbn, 0);

% run solver
s.solve;

% visualize solution
s.gui;

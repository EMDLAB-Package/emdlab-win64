%{
note: calculation of the stator field of a 3-phase induction motor using
emdlab magnetic static solver by exciting stator winding with balanced currents
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% dimensions & parameters
gv_ISD = 74;
gv_OSD = 125;
gv_Lstk = 70;
gv_gap = 0.5;
gv_Dsh = 24;
gv_Ns = 36;
gv_Nr = 26;
gv_wst = 3;
gv_dss = 15;
gv_bs0 = 2.4;
gv_hs0 = 0.5;
gv_stta = 35;
gv_wrt = 4;
gv_drs = 10;
gv_br0 = 2;
gv_hr0 = 0.7;
gv_rtta = 35;
gv_Ntc = 85;
gv_Iph = 1.9;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_stta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_tc4(g, gv_Dsh, gv_ISD-2*gv_gap, gv_Nr, gv_wrt, gv_drs, gv_br0, gv_hr0, gv_rtta, 'rotor', 'rbar', 'rap');

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_gap,gv_OSD/2], [3,gv_gap,3], r, 'linear','extrap');
g.setMeshLengthByRadialFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');
m.joinMeshZones('rb','rbar','rap');
m.setMeshZoneColor('rb',0,192,232);

% add materials
m.addMaterial('m530', emdlab_mlib_es_M530_50A);
m.addMaterial('copper', emdlab_mlib_copper);
m.addMaterial('aluminium', emdlab_mlib_aluminium);

% set materials
m.setMaterial('rotor','m530');
m.setMaterial('stator','m530');
m.setMaterial('sc','copper');
m.setMaterial('rb','aluminium');

% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns)
m.aux_cmxjcrj('sap',gv_Ns)
m.aux_cmxjcr('sc',gv_Ns)
m.aux_cmxjcrj('rotor',gv_Nr)
m.aux_cmxjcr('rb',gv_Nr)

% add circular air gap
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_gap,0,0,gv_ISD/2,2);

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% define windings
s.defineCoil('phaseA');
s.defineCoil('phaseB');
s.defineCoil('phaseC');
pA = [1,2,3,10,11,12]; pA = [pA, pA+18];
ntc = gv_Ntc*[-1,-1,-1,1,1,1]; ntc = [ntc,ntc];
pB = pA + 6;
pC = pA + 3;
s.addMeshZones2Coil('phaseA', 'sc' + string(pA), ntc, 0.15);
s.addMeshZones2Coil('phaseB', 'sc' + string(pB), ntc, 0.15);
s.addMeshZones2Coil('phaseC', 'sc' + string(pC), -ntc, 0.15);

% set phase currents
s.setCoilCurrent('phaseA', gv_Iph*1.41);
s.setCoilCurrent('phaseB', -gv_Iph*1.41/2);
s.setCoilCurrent('phaseC', -gv_Iph*1.41/2);

% apply boundary conditions: zero vector potential on all boundary nodes
s.setAzBC(m.getfbn, 0);

% run solver
s.solve;

% visualize solution
s.gui;
s.plotBrBtOnCircle(0, 0, gv_ISD/2-gv_gap/2, 1000);

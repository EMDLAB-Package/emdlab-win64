%{
note: calculation of the stator field of a three-phase synchronous reluctance motor
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 75;
gv_OSD = 125;
gv_Lstk = 70;
gv_Ns = 36;
gv_p = 4;
gv_wst = 3;
gv_dss = 15;
gv_bs0 = 2.4;
gv_hs0 = 0.6;
gv_tta = 25;
gv_Dsh = 24;
gv_Ntc = 85;
gv_g = 0.4;
gv_Iph = 1.9;
gv_Nfb = 4;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_synrm1(g, gv_Dsh, gv_ISD-2*gv_g, gv_p, 0.45, ones(1,gv_Nfb)*0.8,  ones(1,gv_Nfb-1)*0.8, 0.6, 'rotor', 'fb')

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [2,0.4,3], r, 'linear','extrap');
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
m.aux_unify('fb')
m.aux_cmxjcrj('fb',gv_p)
m.aux_cmxjcr('sc',gv_Ns)

% rotate mesh zones
for mz = ["rotor", "fb"]
    m.rotateMeshZone(mz,pi/4-pi/gv_Ns);
end

% generate air gap mesh
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,2);

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

% apply boundary conditions
s.setAzBC(m.getfbn, 0);

% run solver
s.solve(100);

% visualize solution
s.gui;

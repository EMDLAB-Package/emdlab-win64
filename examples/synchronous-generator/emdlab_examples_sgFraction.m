%{
note: calculation of the rotor field of a synchronous generator
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% dimensions & parameters
gv_Dsh = 60; % shaft diameter [mm]
gv_ISD = 170; % inner stator diamater [mm]
gv_OSD = 250; % outer stator diameter [mm]
gv_Lstk = 150; % stack length [mm]
gv_gap = 1; % airgap thickness [mm]
gv_gg = 2.5;
gv_Ns = 48; % number of stator slots
gv_p = 6; % number of rotor poles
gv_wst = 6; % width of stator tooth [mm]
gv_dss = 25; % depth of stator slot [mm]
gv_bs0 = 2; % stator slot opening [mm] 
gv_hs0 = 1.5; % stator slot opening height [mm]
gv_tta = 15; % tooth tip angle of stator slot [deg]

gv_wpole = 30; % magnet depth [mm]
gv_Nagl = 2;
gv_embrace = 0.75; % magnet embrace ratio -> lower than one
gv_Hc = -922100; % magnet coercive force [A/m]

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_tc_wfsm1(g, gv_Dsh, gv_ISD, gv_gap, gv_gg, gv_Nagl, gv_p, gv_embrace, gv_wpole, 5, 0.7, 0.9, 'rotor', 'fcoil', 'rap');

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [3,0.5,2], r, 'linear','extrap');
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
m.setMaterial('fcoil','copper');

% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns/6,2*pi/gv_Ns)
m.aux_cmxjcrj('sap',gv_Ns/6,2*pi/gv_Ns)
m.aux_cmxjcrj('rotor',gv_p/6,2*pi/gv_p)
m.aux_unify('rap');
m.aux_cmxjcrj('rap',gv_p/6,2*pi/gv_p)
m.aux_cmxcr('sc',gv_Ns/6,2*pi/gv_Ns)
m.aux_cmxcr('fcoil',gv_p/6,2*pi/gv_p)

% rotate mesh zones
m.rotateMeshZones(["stator", "sap", "sc1"+string(1:8), "sc2"+string(1:8)], pi/gv_Ns);
m.rotateMeshZones(["rotor", "fcoil11", "fcoil21", "rap"], pi/gv_p);

% generate air gap mesh
m.aux_addArcAirGap('ag',0,0,gv_ISD/2-(gv_Nagl-1)*gv_gap/gv_Nagl,gv_ISD/2,gv_Nagl-1);

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

s.defineCoil('fcoil');
for i = 1:gv_p/6
    if rem(i,2), dir = [1,-1]; else, dir = [-1,1]; end
    s.addMeshZones2Coil('fcoil', 'fcoil1' + string(i), dir(1));
    s.addMeshZones2Coil('fcoil', 'fcoil2' + string(i), dir(2));
end

% set phase currents
s.setCoilCurrent('fcoil', 600); 

% apply boundary conditions
kr = m.getfbnioc([0,0],gv_OSD/2);
ks = m.getfbnioc([0,0],gv_Dsh/2);
index = [ks;kr];
s.setAzBC(index,0);
fb_index = setdiff(m.getfbn,index);
[km,ks] = m.splitPeriodic(fb_index, 2*pi/gv_p);
s.setOddPeriodicBC(km,ks);

% run solver
s.solve();

% visualize solution
m.showg;
s.plotBmagF(16,'rotor','stator');
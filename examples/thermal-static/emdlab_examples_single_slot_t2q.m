% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 73; % inner stator diameter [mm]
gv_OSD = 125; % outer stator diameter [mm]
gv_Lstk = 70; % stack length [mm]
gv_Ns = 36; % number of stator slots
gv_wst = 3.3; % width of stator tooth
gv_dss = 18; % depth of stator slot [mm]
gv_th = 4; % thickness of housing [mm]
meshSize = 1; % mesh size [mm]

% dependents
gv_alpha_s = 2*pi/gv_Ns; % stator slot pitch angle [rad]
gamma_1 = gv_alpha_s/2 - asin(gv_wst/gv_ISD);
gamma_2 = gv_alpha_s/2 - asin(gv_wst/(gv_ISD+2*gv_dss));

% get 2d geometry data base
g = emdlab_g2d_db;

% add points
g.addPoint(gv_ISD/2,0);
g.addPoint(gv_ISD/2+gv_dss,0);
g.addPoint(gv_OSD/2,0);
g.addPoint(gv_OSD/2+gv_th,0);
g.addPoint((gv_ISD/2)*cos(gamma_1),(gv_ISD/2)*sin(gamma_1));
g.addPoint((gv_ISD/2+gv_dss)*cos(gamma_2),(gv_ISD/2+gv_dss)*sin(gamma_2));
g.addPoint((gv_OSD/2)*cos(gamma_2),(gv_OSD/2)*sin(gamma_2));
g.addPoint((gv_OSD/2+gv_th)*cos(gamma_2),(gv_OSD/2+gv_th)*sin(gamma_2));
g.addPoint((gv_ISD/2)*cos(gv_alpha_s/2),(gv_ISD/2)*sin(gv_alpha_s/2));
g.addPoint((gv_ISD/2+gv_dss)*cos(gv_alpha_s/2),(gv_ISD/2+gv_dss)*sin(gv_alpha_s/2));
g.addPoint((gv_OSD/2)*cos(gv_alpha_s/2),(gv_OSD/2)*sin(gv_alpha_s/2));
g.addPoint((gv_OSD/2+gv_th)*cos(gv_alpha_s/2),(gv_OSD/2+gv_th)*sin(gv_alpha_s/2));
g.addPoint(0,0);

% add edges
g.addSegment(1,2);
g.addSegment(2,3);
g.addSegment(3,4);
g.addSegment(5,6);
g.addSegment(6,7);
g.addSegment(7,8);
g.addSegment(9,10);
g.addSegment(10,11);
g.addSegment(11,12);
g.addArc(13,1,5,1);
g.addArc(13,5,9,1);
g.addArc(13,2,6,1);
g.addArc(13,6,10,1);
g.addArc(13,3,7,1);
g.addArc(13,7,11,1);
g.addArc(13,4,8,1);
g.addArc(13,8,12,1);

g.addFace('copper', g.addLoop(1,12,-4,-10));
g.addFace('sz1', g.addLoop(4,13,-7,-11));
g.addFace('sz2', g.addLoop(2,14,-5,-12));
g.addFace('sz3', g.addLoop(5,15,-8,-13));
g.addFace('hz1', g.addLoop(3,16,-6,-14));
g.addFace('hz2', g.addLoop(6,17,-9,-15));

g.setMeshMaxLength(meshSize);
m = g.generateMesh('gmsh');
m.joinMeshZones('stator', 'sz1', 'sz2',' sz3');
m.aux_cmxjcrj('stator',1)
m.joinMeshZones('housing', 'hz1', 'hz2');
m.aux_cmxjcrj('housing',1)
m.aux_cmxjcrj('copper',1)
m.showce
% get mesh database
m = m.aux_c2qm;
m.setMeshZoneColor('stator', 190, 190, 190);
m.setMeshZoneColor('housing',163,73,164);
m.setMeshZoneColor('copper',255,137,39);

% add materials
m.addMaterial('copper', emdlab_mlib_copper);
m.addMaterial('iron', emdlab_mlib_iron);
m.addMaterial('aluminium', emdlab_mlib_aluminium);

% set materials
m.mts.iron.setThermalConductivity(50);
m.setMaterial('stator', 'iron');
m.mts.copper.setThermalConductivity(10);
m.setMaterial('copper', 'copper');
m.mts.aluminium.setThermalConductivity(200);
m.setMaterial('housing', 'aluminium');

% add solver
s = emdlab_solvers_ts2d_tn_qm(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set boundary conditions
idx = m.getEdgeIndicesInCircle(0,0,gv_ISD/2);
s.addHeatFluxBC('inner_surface', idx, 500);

idx = m.getEdgeIndicesOnCircle(0,0,gv_OSD/2 + gv_th,(gv_OSD/2 + gv_th)*(1-cos(asin(meshSize/(gv_OSD/2 + gv_th)))));
s.addConvectionBC('outer_surface', idx, 10, 25);

% solve and plot results
s.solve;
m.showmzs;
s.plotAverageTemperature(20);
s.plotTemperature(20);
fprintf('Tmin = %.4f\n', min(s.results.T));
fprintf('Tmax = %.4f\n', max(s.results.T));
fprintf('Tavg = %.4f\n', s.getAverageTemperature);

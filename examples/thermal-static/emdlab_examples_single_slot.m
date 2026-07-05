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

Nr1 = ceil(g.getEdgeLength(1)/meshSize);
Nr2 = ceil(g.getEdgeLength(2)/meshSize);
Nr3 = ceil(g.getEdgeLength(3)/meshSize);
Nt1 = ceil(g.getEdgeLength(16)/meshSize);
Nt2 = ceil(g.getEdgeLength(17)/meshSize);

% get mesh database
m = emdlab_m2d_qmdb;
m.addMeshZone('copper', g.getQMeshByEdges(1,12,-4,-10,Nr1,Nt1));
m.addMeshZone('sz1', g.getQMeshByEdges(4,13,-7,-11,Nr1,Nt2));
m.addMeshZone('sz2', g.getQMeshByEdges(2,14,-5,-12,Nr2,Nt1));
m.addMeshZone('sz3', g.getQMeshByEdges(5,15,-8,-13,Nr2,Nt2));
m.addMeshZone('hz1', g.getQMeshByEdges(3,16,-6,-14,Nr3,Nt1));
m.addMeshZone('hz2', g.getQMeshByEdges(6,17,-9,-15,Nr3,Nt2));

m.joinMeshZones('stator', 'sz1', 'sz2',' sz3');
m.aux_cmxjcrj('stator',1)
m.setMeshZoneColor('stator', 190, 190, 190);

m.joinMeshZones('housing', 'hz1', 'hz2');
m.aux_cmxjcrj('housing',1)
m.setMeshZoneColor('housing',163,73,164);

m.aux_cmxjcrj('copper',1)
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
s = emdlab_solvers_ts2d_lptn_qm(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set boundary conditions
idx = m.getEdgeIndicesOnCircle(0,0,gv_ISD/2);
s.addHeatFluxBC('inner_surface', idx, 500);

idx = m.getEdgeIndicesOnCircle(0,0,gv_OSD/2 + gv_th);
s.addConvectionBC('outer_surface', idx, 10, 25);

% solve and plot results
s.solve;
m.showmzs;
s.plotAverageTemperature(20);
mean(s.results.T)
s.plotTemperature(20);

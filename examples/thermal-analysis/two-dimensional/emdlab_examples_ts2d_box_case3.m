%{
Solving 2D heat diffusion equation in box with
fixed temperature at left edge and convection for 
the rest
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% inputs
W = 1; % width of the box
H = 1; % height of the box
Z = 1; % depth of the problem
meshSize = 1/20; % maximum mesh size

% define geometry
g = emdlab_g2d_db;
g.addRectangleLoop(0,0,W,H);

% construct quadrilateral mesh
m = emdlab_m2d_qmdb();
m.addMeshZone('z1', g.getQMeshByEdges(1,2,3,4,ceil(W/meshSize),ceil(H/meshSize)));

% add materials
m.addMaterial('copper', emdlab_mlib_copper);
m.setMaterial('z1', 'copper');
m.setMeshZoneColor('z1', 0, 255, 255);
m.mts.copper.setThermalConductivity(1);

% add solver
s = emdlab_solvers_ts2d_tn(m);
s.setLengthUnit('m');
s.setDepth(1);

% set boundary conditions
idx_left = m.getEdgeIndicesOnLineP0U(0,0,0,1);
s.addFixedTemperatureBC('left', idx_left, 1);

idx_bottom = m.getEdgeIndicesOnLineP0U(0,0,1,0);
idx_top = m.getEdgeIndicesOnLineP0U(0,H,1,0);
idx_right = m.getEdgeIndicesOnLineP0U(W,0,0,1);
s.addConvectionBC('rest', [idx_bottom;idx_top;idx_right], 10, 0);

% solve and plot results
s.solve;
s.plotThermalNetwork;
s.plotTemperature(10);
fprintf('Tmin = %.4f\n', min(s.results.T));
fprintf('Tmax = %.4f\n', max(s.results.T));
fprintf('Tavg = %.4f\n', s.getAverageTemperature);

fprintf('q_left = %.4f\n', s.calculateNetHeatCrossingBoundaryEdges(idx_left));
fprintf('q_right = %.4f\n', s.calculateNetHeatCrossingBoundaryEdges(idx_right));
fprintf('q_top = %.4f\n', s.calculateNetHeatCrossingBoundaryEdges(idx_top));
fprintf('q_bottom = %.4f\n', s.calculateNetHeatCrossingBoundaryEdges(idx_bottom));
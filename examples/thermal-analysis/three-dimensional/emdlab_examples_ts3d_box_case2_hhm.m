%{
Solving 3D heat diffusion equation in box with
input heat-flux at left face and zero temperature for 
the rest -> using hexahedron mesh
Tavg = 3.43612
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
meshSize = 0.05; % maximum mesh size

% define geometry
g = emdlab_g2d_db;
g.addRectangleLoop(0,0,W,H);

% construct quadrilateral mesh
qm = emdlab_m2d_qmdb();
qm.addMeshZone('z1', g.getQMeshByEdges(1,2,3,4,ceil(W/meshSize),ceil(H/meshSize)));

% extrude quadrilateral mesh to generate hexahedron mesh
m = emdlab_m3d_hhmdb;
m.addMeshZone('z1', qm.mzs.z1.getExtrude(linspace(0,1,ceil(1/meshSize))));

% add & set materials
m.addMaterial('copper', emdlab_mlib_copper);
m.setMaterial('z1', 'copper');
m.mts.copper.setThermalConductivity([1,1,1]);

% define the solver
s = emdlab_solvers_ts3d_tn(m);

% set left face boundary condition
left_idx = m.getFacetIndicesOnPlane([0,0,0],[1,0,0]);
s.addHeatFluxBC('left', left_idx, 100);

% set boundary condition for the rest faces
rest_idx = setdiff(m.getfbf, left_idx);
s.addFixedTemperatureBC('rest', rest_idx, 0);

% solve & plot results
s.solve
s.plotTemperature;
fprintf('Tmin = %.4f\n', s.getMinimumTemperature);
fprintf('Tmax = %.4f\n', s.getMaximumTemperature);
fprintf('Tavg = %.4f\n', s.getAverageTemperature);
fprintf('Qin = %.4f\n', s.calculateNetHeatCrossingBoundaryFacets(left_idx));
fprintf('Qout = %.4f\n', s.calculateNetHeatCrossingBoundaryFacets(rest_idx));
right_idx = m.getFacetIndicesOnPlane([1,0,0],[1,0,0]);
fprintf('Qright_face = %.4f\n', s.calculateNetHeatCrossingBoundaryFacets(right_idx));

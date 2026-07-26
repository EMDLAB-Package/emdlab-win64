%{
Solving 3D heat diffusion equation in box with
non-zero temperature at left face and zero temperature for 
the rest -> using hexahedron mesh
Tavg = 0.8905
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
s.addFixedTemperatureBC('left', left_idx, @(x,y,z) 10*sin(pi*y)*sin(pi*z));

% set boundary condition for the rest faces
rest_idx = setdiff(m.getfbf, left_idx);
s.addFixedTemperatureBC('rest', rest_idx, 0);

% solve & plot results
s.solve
s.plotAverageTemperature;
fprintf('Tmin = %.4f\n', min(s.results.T));
fprintf('Tmax = %.4f\n', max(s.results.T));
fprintf('Tavg = %.4f\n', s.getAverageTemperature);

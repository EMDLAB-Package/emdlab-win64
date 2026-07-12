%{
Solving 2D heat diffusion equation in circle with
zero temperature for all boundaries with
internal heat generation
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% inputs
R = 1; % circle radius
Z = 1; % depth of the problem
meshSize = 0.2; % maximum mesh size

% % construct quadrilateral mesh
m = emdlab_m2d_qmdb();
m.addMeshZone('z1', emdlab_m2d_gqm4circle(0,0,R,meshSize));

% add materials
m.addMaterial('copper', emdlab_mlib_copper);
m.setMaterial('z1', 'copper');
m.setMeshZoneColor('z1', 0, 255, 255);
m.mts.copper.setThermalConductivity(1);

% add solver
s = emdlab_solvers_ts2d_tn_qm(m);
s.setLengthUnit('m');
s.setDepth(Z);

% set boundary condition & internal heat generation
s.addFixedTemperatureBC('zeroT', m.getfbe, 0);
s.addInternalHeatSource('hg1', 'z1', 1000);

% solve and plot results
s.solve;
s.plotThermalNetwork;
s.plotTemperature(10);
fprintf('Tmin = %.4f\n', min(s.results.T));
fprintf('Tmax = %.4f\n', max(s.results.T));
fprintf('Tavg = %.4f\n', s.getAverageTemperature);
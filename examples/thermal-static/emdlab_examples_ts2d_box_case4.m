%{
Solving 2D heat diffusion equation in box with
zero temperature for all boundaries with
internal heat generation
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
meshSize = 1/100; % maximum mesh size

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
m.showmzs;
m.mts.copper.ThermalConductivity.value = 1;

% add solver
s = emdlab_solvers_ts2d_lptn_qm(m);
s.setLengthUnit('m');
s.setDepth(Z);

% set boundary condition & internal heat generation
s.addFixedTemperatureBC('left', m.getfbe, 0);
s.addInternalHeatSource('hg1', 'z1', 1000);

% solve and plot results
s.solve;
s.plotAverageTemperature(20);
mean(s.results.T)

%{
note: calculation of the rotor field of a single-rotor single-stator
axial-flux permanent magnet motor
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 150;
gv_OSD = 250;
gv_Ns = 24;
gv_p = 4;
gv_wss = 12;
gv_dss = 40;
gv_dsy = 10;
gv_Ntc = 85;
gv_gap = 1;
gv_Iph = 1.9;
gv_dm = 6;
gv_dry = 10;
gv_embrace = 0.8;
gv_Hc = -922100;

% dependent variables
gv_Lstk = (gv_OSD-gv_ISD)/2;
gv_D = (gv_OSD+gv_ISD)/2;
gv_L = pi*gv_D;
gv_taus = gv_L/gv_Ns;
gv_wst = gv_taus - gv_wss;
gv_taup = gv_L/gv_p;

% define geometry data base
g = emdlab_g2d_db;

% add geometries from templates
emdlab_g2d_lib_tc110(g,gv_wst,gv_wss,gv_dss,5,2,25,0,0,0,'stator_1','sc','sap');
% emdlab_g2d_lib_rm_spm50(g,gv_taup,gv_dry,gv_dm,0,gv_dss+gv_dsy+gv_gap+gv_dm+gv_dry,1,'rotor','magnet');

% setting the wireframe mesh by mesh size function
f_mesh = @(y) interp1([0,gv_dsy+gv_dss,gv_gap+gv_dm+gv_dry+gv_dsy+gv_dss],[gv_dsy/3,gv_gap/2,gv_dry/3],y, 'linear','extrap');
g.setMeshLengthByYFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');


% m.aux_cmyjcshj('rotor',1,gv_taup,0)
% m.aux_cmyjcshj('magnet',1,gv_taup,0)
m.copyMirrorMeshZone('stator_2','stator_1',[1,0]);
m.joinMeshZones('stator','stator_2','stator_1')
m.aux_cmyjcshj('stator',6,gv_taus,0)
m.aux_cmyjcshj('sc',6,gv_taus,0)
m.aux_cmyjcshj('sap',6,gv_taus,0)
m.shiftMeshZone('stator',-gv_taup/2+gv_taus/2,0);
m.shiftMeshZone('sc',-gv_taup/2+gv_taus/2,0);
m.shiftMeshZone('sap',-gv_taup/2+gv_taus/2,0);
m.showg
% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design parameters
gv_ISD = 100;
gv_OSD = 150;
gv_g = 0.7;
gv_Lstk = 100;
gv_Dsh = 55;
gv_Ns = 48;
gv_Nc = 6;
gv_p = 8;
gv_Hc = -922100;

% dependent variables
alpha_p = 2*pi/gv_p;
alpha_s = 2*pi/gv_Ns;

% generator geometry
g = emdlab_g2d_db;
% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [2,0.5,2], r, 'linear', 'extrap');
emdlab_g2d_lib_tc9(g,gv_ISD,gv_OSD,gv_Ns,gv_Nc,3,1.8,1.5,0.8,15,0.3,0.3,'stator','sc','sap');
emdlab_g2d_lib_rm_ipm13(g, gv_Dsh, gv_ISD-2*gv_g, gv_p, 0.25, [0.8,1], [1]*.8, [0.6,0.7], [1,1]*0.5, 0.6, 'rotor', 'magnet', 'rap');
g.setMeshLengthByRadialFunction(f_mesh)
m = g.generateMesh('mg0');

% add materials
m.addMaterial('m330', emdlab_mlib_es_M300_35A);
m.addMaterial('copper', emdlab_mlib_copper);

% set materials
m.setMaterial('rotor','m330');
m.setMaterial('stator','m330');

% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns/gv_p, alpha_s);
m.aux_cmxjcrj('sap',gv_Ns/gv_p, alpha_s);
for i = 1:gv_Nc
    m.aux_cmxjcrj('sc'+string(i),gv_Ns/gv_p, alpha_s);
end
m.aux_cmxjcrj('rotor',1, alpha_p);
m.aux_unify('rap');
m.aux_cmxjcrj('rap',1, alpha_p);
m.joinMeshZones('magnetc', 'magnetc'+string(1:2))
m.aux_cmxjcr('magnetc',1, alpha_p);
m.joinMeshZones('magnets', 'magnets'+string(1:2))
m.aux_cmxcr('magnets',1, alpha_p);
m.aux_addArcAirGap('ag',0,0, gv_ISD/2-gv_g,gv_ISD/2,2);

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set magnetization
u1 = [cos(alpha_p/2),sin(alpha_p/2)];
u1s = emdlab_g2d_rotatePoints(u1,-pi/2);
u2s = emdlab_g2d_rotatePoints([u1(1),-u1(2)],pi/2);
s.setMagnetization('magnetc1',gv_Hc,[1,0]);
s.setMagnetization('magnets11',gv_Hc,u1s);
s.setMagnetization('magnets21',gv_Hc,u2s);

% apply boundary conditions
kr = m.getfbnioc([0,0],gv_OSD/2);
ks = m.getfbnioc([0,0],gv_Dsh/2);
index = [ks;kr];
s.setAzBC(index,0);
fb_index = setdiff(m.getfbn,index);
[km,ks] = m.splitPeriodic(fb_index, 2*pi/gv_p);
s.setOddPeriodicBC(km',ks');

% solve and plot results
s.assignEdata(20)
s.solve;
s.gui;
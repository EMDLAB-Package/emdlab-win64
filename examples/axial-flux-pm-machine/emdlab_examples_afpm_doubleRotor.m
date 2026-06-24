%{
note: calculation of the rotor field of a double-rotor
axial-flux permanent magnet motor
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 120;
gv_OSD = 160;
gv_Ns = 12;
gv_p = 10;
gv_wss = 16;
gv_dss = 40;
gv_Ntc = 85;
gv_gap = 1;
gv_Iph = 1.9;
gv_dm = 5;
gv_dry = 10;
gv_embrace = 0.8;
gv_Laxial = gv_dss + 2*(gv_gap+gv_dm+gv_dry);
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
emdlab_g2d_lib_tc110(g,gv_wst,gv_wss,gv_dss,5,2,25,0,0,0,'s1','sc1','sap1');
emdlab_g2d_lib_rm_spm102(g,gv_taup,gv_dry,gv_dm,0.85,0,gv_dss/2+gv_gap+gv_dm+gv_dry,1 ,'rotorFront','magnetFront','rapFront');

% setting the wireframe mesh by mesh size function
f_mesh = @(y) interp1([0,gv_dss/2,gv_dss/2+gv_gap+gv_dry+gv_dm],[gv_wst/4,gv_gap,gv_dry/3],y, 'linear','extrap');
g.setMeshLengthByYFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');

% add materials
m.addMaterial('m530', emdlab_mlib_es_M530_50A);
m.addMaterial('copper', emdlab_mlib_copper);

% set mesh zone materials
m.setMaterial('s1','m530');
m.setMaterial('rotorFront','m530');
m.setMaterial('sc1','copper');

m.copyMirrorMeshZone('rotorRear','rotorFront',[1,0]);
m.copyMirrorMeshZone('magnetRear','magnetFront',[1,0]);
m.copyMirrorMeshZone('rapRear','rapFront',[1,0]);

% generate full mesh
m.aux_cmyjcshj('rotorRear',gv_p,gv_taup,0)
m.aux_cmyjcshj('rapRear',gv_p,gv_taup,0)
m.aux_cmyjcsh('magnetRear',gv_p,gv_taup,0)
m.aux_cmyjcshj('rotorFront',gv_p,gv_taup,0)
m.aux_cmyjcshj('rapFront',gv_p,gv_taup,0)
m.aux_cmyjcsh('magnetFront',gv_p,gv_taup,0)

m.copyMirrorMeshZone('s2','s1',[1,0]);
m.joinMeshZones('stator','s1','s2');
m.aux_cmyjcshj('stator',gv_Ns,gv_taus,0)
m.copyMirrorMeshZone('sc2','sc1',[1,0]);
m.joinMeshZones('sca','sc1','sc2');
m.aux_cmycsh('sca',gv_Ns,gv_taus,0)
m.copyMirrorMeshZone('sap2','sap1',[1,0]);
m.joinMeshZones('sap','sap1','sap2');
m.aux_cmyjcshj('sap',gv_Ns,gv_taus,0);

% add airgap mesh
m.aux_addLineAirGap('ag1',0,gv_dss/2+gv_gap/2,1,0,gv_gap,2)
m.aux_addLineAirGap('ag2',0,-gv_dss/2-gv_gap/2,1,0,gv_gap,2)

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set magnetization
for i = 1:2:gv_p
s.setMagnetization(['magnetFront',num2str(i)],gv_Hc,[0,1]);
m.setMeshZoneColor(['magnetFront',num2str(i)],0,125,223);
s.setMagnetization(['magnetRear',num2str(i)],gv_Hc,[0,1]);
m.setMeshZoneColor(['magnetRear',num2str(i)],0,125,223);
end
for i = 2:2:gv_p
s.setMagnetization(['magnetFront',num2str(i)],gv_Hc,[0,-1]);
m.setMeshZoneColor(['magnetFront',num2str(i)],255,70,70);
s.setMagnetization(['magnetRear',num2str(i)],gv_Hc,[0,-1]);
m.setMeshZoneColor(['magnetRear',num2str(i)],255,70,70);
end

% apply boundary conditions
k1 = m.getfbniol_p0x([0,gv_Laxial/2]);
k2 = m.getfbniol_p0x([0,-gv_Laxial/2]);
s.setAzBC([k1;k2],0);
[km,ks] = m.splitShift(setdiff(m.getfbn,[k1;k2]), [gv_L,0]);
s.setEvenPeriodicBC(km,ks);

% run solver
s.solve;
s.plotBmagF(8, 'stator', 'rotorRear', 'rotorFront');
cb = findall(gcf,'Type','ColorBar');
cb.Location = 'southoutside';
zoom(2);
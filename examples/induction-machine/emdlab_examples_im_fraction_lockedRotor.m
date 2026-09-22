%{
note: transient voltage-fed locked rotor simulation of an induction motor -> fraction model
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% dimensions & parameters
gv_ISD = 74;
gv_OSD = 125;
gv_Lstk = 70;
gv_gap = 0.5;
gv_Dsh = 24;
gv_Ns = 36;
gv_Nr = 26;
gv_wst = 3;
gv_dss = 15;
gv_bs0 = 2.4;
gv_hs0 = 0.5;
gv_stta = 35;
gv_wrt = 4;
gv_drs = 10;
gv_br0 = 2;
gv_hr0 = 0.7;
gv_rtta = 35;
gv_Ntc = 85;
gv_Iph = 1.9;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_stta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_tc4(g, gv_Dsh, gv_ISD-2*gv_gap, gv_Nr, gv_wrt, gv_drs, gv_br0, gv_hr0, gv_rtta, 'rotor', 'rbar', 'rap');

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_gap,gv_OSD/2], [3,gv_gap,3], r, 'linear','extrap');
g.setMeshLengthByRadialFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');
m.joinMeshZones('rb','rbar','rap');
m.setMeshZoneColor('rb',0,192,232);

% add materials
m.addMaterial('m530', emdlab_mlib_es_M530_50A);
m.addMaterial('copper', emdlab_mlib_copper);
m.addMaterial('aluminium', emdlab_mlib_aluminium);

% set materials
m.setMaterial('rotor','m530');
m.setMaterial('stator','m530');
m.setMaterial('sc','copper');
m.setMaterial('rb','aluminium');

% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns/2,2*pi/gv_Ns)
m.aux_cmxjcrj('sap',gv_Ns/2,2*pi/gv_Ns)
m.aux_cmxjcr('sc',gv_Ns/2,2*pi/gv_Ns)
m.aux_cmxjcrj('rotor',gv_Nr/2,2*pi/gv_Nr)
m.aux_cmxjcr('rb',gv_Nr/2,2*pi/gv_Nr)

% add circular air gap
m.aux_addArcAirGap('ag',0,0,gv_ISD/2-gv_gap,gv_ISD/2,2);

% getting an instance of solver object
s = emdlab_solvers_mt2d_tl3_ihnlwtm(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% define windings
s.defineCoil('phaseA','voltage','stranded');
s.defineCoil('phaseB','voltage','stranded');
s.defineCoil('phaseC','voltage','stranded');
% index of each coil assocciating with each winding
pA = [1,2,3,10,11,12];
ntc = gv_Ntc*[-1,-1,-1,1,1,1];
pB = pA + 6;
pC = pA + 3;

% assingation of mesh zones to windings
s.addMeshZones2Coil('phaseA', 'sc' + string(pA), ntc, 0.15);
s.addMeshZones2Coil('phaseB', 'sc' + string(pB), ntc, 0.15);
s.addMeshZones2Coil('phaseC', 'sc' + string(pC), -ntc, 0.15);

% set phase voltages
s.setCoilVoltage('phaseA', @(t) 0.5*325*sin(2*pi*50*t));
s.setCoilVoltage('phaseB', @(t) 0.5*325*sin(2*pi*50*t - 2*pi/3));
s.setCoilVoltage('phaseC', @(t) 0.5*325*sin(2*pi*50*t + 2*pi/3));

% apply boundary conditions
kr = m.getnIndexOnCircle([0,0],gv_OSD/2);
ks = m.getnIndexOnCircle([0,0],gv_Dsh/2);
index = [ks;kr];
s.setAzBC(index,0);
fb_index = setdiff(m.getfbn,index);
[km,ks] = m.splitPeriodic(fb_index, pi);
s.setEvenPeriodicBC(km,ks);

s.defineStarConnection('y1', 'phaseA', 'phaseB', 'phaseC')
s.defineCage('cage_1', 'rb' + string(1:13));

% solve for 250ms, time step is 1ms
% a vector to store the value of torque
te = zeros(1,250);

% run solver step by step
for i = 1:250
    s.solveForOneTimeStep(1e-3);
    te(i) = s.evalTorqueByArkkio('ag', gv_gap);
end

% plot results
g.showSketch;
m.showmzs;
s.plotCoilCurrents
s.plotCoilFluxLinkages
figure;
plot(te);
xlabel('Time [ms]');
ylabel('Elecrtic Torque [Nm]');
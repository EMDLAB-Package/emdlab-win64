%{
note: back-emf simulation of a surface-mounted PMSM
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 74;
gv_OSD = 125;
gv_Lstk = 70;
gv_Ns = 36;
gv_p = 4;
gv_wst = 3;
gv_dss = 15;
gv_bs0 = 2.4;
gv_hs0 = 0.6;
gv_tta = 25;
gv_Dsh = 28;
gv_dm = 3;
gv_g = 1;
gv_Ntc = 65;
gv_embrace = 0.82;
gv_Hc = -922100;

% define geometry data base
g = emdlab_g2d_db;

% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_dm-gv_g,gv_ISD/2,gv_OSD/2], [4,1,gv_g,3], r, 'linear','extrap');

% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_spm1(g, gv_Dsh, gv_ISD-2*gv_g, gv_p, gv_dm, gv_embrace, 'rotor', 'magnet', 'rap')

% set mesh density function
g.setMeshLengthByRadialFunction(f_mesh);
m = g.generateMesh('mg0');

% add materials
m.addMaterial('m330', emdlab_mlib_es_M330_35A);
m.addMaterial('copper', emdlab_mlib_copper);

% set materials
m.setMaterial('rotor','m330');
m.setMaterial('stator','m330');
m.setMaterial('sc','copper');

% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns)
m.aux_cmxjcrj('sap',gv_Ns)
m.aux_cmxjcrj('rotor',gv_p)
m.aux_cmxjcrj('rap',gv_p)
m.aux_cmxjcr('sc',gv_Ns)
m.aux_cmxjcr('magnet',gv_p)

% generate air gap mesh
m.aux_addCircularAirGapInterface('ag',0,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,3, 'inner')

% getting an instance of solver object
s = emdlab_solvers_mt2d_tl3_ihnlwm(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set magnetization
for i = 1:2:gv_p
    s.setMagnetization(['magnet',num2str(i)],gv_Hc,'r');
    m.setMeshZoneColor(['magnet',num2str(i)],0,125,223);
end
for i = 2:2:gv_p
    s.setMagnetization(['magnet',num2str(i)],gv_Hc,'-r');
    m.setMeshZoneColor(['magnet',num2str(i)],255,70,70);
end

% define windings
s.defineCoil('phaseA','current','stranded');
s.defineCoil('phaseB','current','stranded');
s.defineCoil('phaseC','current','stranded');
% index of each coil assocciating with each winding
pAn = [1,2,3]; pAn = [pAn, pAn+18];
pAp = [10,11,12]; pAp = [pAp, pAp+18];
pBp = pAp + 6;
pBn = pAn + 6;
pCp = pAn + 3;
pCn = pAp + 3;
% assingation of mesh zones to windings
for i = 1:6
    s.addMeshZone2Coil('phaseA', ['sc',num2str(pAp(i))], gv_Ntc, 1, 0.15);
    s.addMeshZone2Coil('phaseA', ['sc',num2str(pAn(i))], gv_Ntc, -1, 0.15);
    s.addMeshZone2Coil('phaseB', ['sc',num2str(pBp(i))], gv_Ntc, 1, 0.15);
    s.addMeshZone2Coil('phaseB', ['sc',num2str(pBn(i))], gv_Ntc, -1, 0.15);
    s.addMeshZone2Coil('phaseC', ['sc',num2str(pCp(i))], gv_Ntc, 1, 0.15);
    s.addMeshZone2Coil('phaseC', ['sc',num2str(pCn(i))], gv_Ntc, -1, 0.15);
end
% define moving region
s.defineMovingRegion('moving_1', ["rotor", "rap", "magnet"+string(1:gv_p)], 'ag');
% apply boundary conditions
s.setAzBC(m.getfbn, 0);
% disable solver monitor
s.setMonitor(0);
% synchronous speed
wm = 1500*pi/30;
% total time steps for 20e-3 simulation
Nt = 72; 
% simulation time step size
timeStep = 20e-3/Nt;
te = zeros(1,Nt);
te_mst = zeros(1,Nt);
s.solveForInitialConditions;
for i = 1:Nt    
    % calculate torque
    te(i) = s.evalTorqueByArkkio('ag', gv_g);
    te_mst(i) = s.evalTorqueByMST3(0, 0, gv_ISD/2-gv_g/2, gv_g, 500) ;
    % rotate moving region
    s.rotateMovingRegion('moving_1', wm*timeStep);
    s.solveForOneTimeStep(timeStep);
end
s.forcePeriodicCoilVoltages;
% plot results
s.plotCoilFluxLinkages;
s.plotCoilInducedVoltages;
% cogging torque
figure;
hold on; box on;
plot(s.simTime(2:end), te);
plot(s.simTime(2:end), te_mst);
xlabel('Time [s]'); ylabel('Electric torque [Nm]');
% set(gca,'ylim',[0,1.2*max(te)]);
legend('Arkkio', 'MST');
hold off;
% finding alpha_d4A
figure;
hold on; box on;
theta_E = 2*wm*s.simTime*180/pi;
theta_Ei = linspace(theta_E(1),theta_E(end),10000);
fla_i = interp1(theta_E,s.coils.phaseA.fluxLinkage,theta_Ei,'spline');
plot(theta_E, s.coils.phaseA.fluxLinkage,'o', 'Color','b', 'LineWidth',2);
plot(theta_Ei, fla_i,'Color','k');
[fla_max,index] = max(fla_i);
plot(theta_Ei([index,index]), [-fla_max,fla_max], 'r', 'LineStyle','--', 'LineWidth',2);
title(sprintf('alpha_{d4A} = %.2f', theta_Ei(index)));
set(gca, 'xlim', [0,360]);
xlabel('Electrical angle [deg]');
ylabel('Flux linkage [Wb]');

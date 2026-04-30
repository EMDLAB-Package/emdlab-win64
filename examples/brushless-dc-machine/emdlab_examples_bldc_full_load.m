%{
note: full-load current-fed fixed-speed simulation of an inner-rotor bldc motor
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
gv_Ns = 18;
gv_p = 20;
gv_wst = 6;
gv_dss = 15;
gv_bs0 = 5;
gv_hs0 = 1.5;
gv_tta = 15;
gv_Dsh = 28;
gv_dm = 3;
gv_gap = 1;
gv_Ntc = 3;
gv_Is = 180;
gv_rpm = 1500;
gv_embrace = 0.89;
gv_Hc = -922100;

% define geometry data base
g = emdlab_g2d_db;

% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_dm-gv_gap,gv_ISD/2,gv_OSD/2], [5,1,0.5,4], r, 'linear','extrap');

% add geometry from library
emdlab_g2d_lib_tc1(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_spm1(g, gv_Dsh, gv_ISD-2*gv_gap, gv_p, gv_dm, gv_embrace, 'rotor', 'magnet', 'rap');

% setting the wireframe mesh by mesh size function
g.setMeshLengthByRadialFunction(f_mesh);
m = g.generateMesh('mm');

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
m.aux_cmcr('sc',[cos(pi/gv_Ns),sin(pi/gv_Ns)],gv_Ns)
m.aux_cmxjcr('magnet',gv_p)

% add circular air gap
m.aux_addCircularAirGapInterface('ag',0,0,gv_ISD/2-gv_gap,0,0,gv_ISD/2,3);

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

% apply boundary conditions
s.setAzBC(m.getfbn, 0);

% define windings
s.defineCoil('phaseA','current','stranded');
s.defineCoil('phaseB','current','stranded');
s.defineCoil('phaseC','current','stranded');
% index of each coil assocciating with each winding
pA = [1    10    -2   -11    -9   -18];
pB = [7    16    -8   -17   -15    -6];
pC = [13     4   -14    -5    -3   -12];
% assingation of mesh zones to windings
for i = 1:length(pA)

    mz1Name = "sc1" + string(abs(pA(i)));
    mz2Name = "sc2" + string(abs(pA(i)));
    if pA(i)>0
        s.addMeshZone2Coil('phaseA', mz1Name, gv_Ntc, 1);
        s.addMeshZone2Coil('phaseA', mz2Name, gv_Ntc, -1);
    else
        s.addMeshZone2Coil('phaseA', mz1Name, gv_Ntc, -1);
        s.addMeshZone2Coil('phaseA', mz2Name, gv_Ntc, 1);
    end

    mz1Name = "sc1" + string(abs(pB(i)));
    mz2Name = "sc2" + string(abs(pB(i)));
    if pB(i)>0
        s.addMeshZone2Coil('phaseB', mz1Name, gv_Ntc, 1);
        s.addMeshZone2Coil('phaseB', mz2Name, gv_Ntc, -1);
    else
        s.addMeshZone2Coil('phaseB', mz1Name, gv_Ntc, -1);
        s.addMeshZone2Coil('phaseB', mz2Name, gv_Ntc, 1);
    end

    mz1Name = "sc1" + string(abs(pC(i)));
    mz2Name = "sc2" + string(abs(pC(i)));
    if pC(i)>0
        s.addMeshZone2Coil('phaseC', mz1Name, gv_Ntc, 1);
        s.addMeshZone2Coil('phaseC', mz2Name, gv_Ntc, -1);
    else
        s.addMeshZone2Coil('phaseC', mz1Name, gv_Ntc, -1);
        s.addMeshZone2Coil('phaseC', mz2Name, gv_Ntc, 1);
    end

end

% define moving region
s.defineMovingRegion('moving_1', ["rotor", "rap", "magnet"+string(1:gv_p)], 'ag');

% apply boundary conditions
s.setAzBC(m.getfbn, 0);

% synchronous speed
wm = gv_rpm*pi/30;

% phase currents in stator reference frame
alpha_Is_rrf = 90 * pi/180;
alpha_d4A = -280 * pi/180;
alpha_Is_srf = @(t) alpha_Is_rrf + alpha_d4A + (gv_p/2)*wm*t;
id_srf = @(t) gv_Is * cos(alpha_Is_srf(t));
iq_srf = @(t) gv_Is * sin(alpha_Is_srf(t));
ia_srf = @(t) id_srf(t);
ib_srf = @(t) -id_srf(t)/2 + sqrt(3)*iq_srf(t)/2;
ic_srf = @(t) -id_srf(t)/2 - sqrt(3)*iq_srf(t)/2;
s.setCoilCurrent('phaseA', ia_srf);
s.setCoilCurrent('phaseB', ib_srf);
s.setCoilCurrent('phaseC', ic_srf);

% total time steps for one electric period simulation
Nt = 50; 

% simulation time step size
stopTime = 1/((gv_p/2)*wm/(2*pi));
timeStep = stopTime/Nt;

% allocate vector to store torque
te = zeros(1,Nt+1);
te_mst = zeros(1,Nt+1);
s.solveForInitialConditions;
te(1) = s.evalTorqueByArkkio('ag', gv_gap);
te_mst(1) = s.evalTorqueByMST(0,0,gv_ISD/2-gv_gap/2,1000);

for i = 1:Nt    

    % rotate moving region
    s.rotateMovingRegion('moving_1', wm*timeStep);
    s.solveForOneTimeStep(timeStep);

    % calculate torque
    te(i+1) = s.evalTorqueByArkkio('ag', gv_gap);
    te_mst(i+1) = s.evalTorqueByMST3(0,0,gv_ISD/2-gv_gap/2,gv_gap,1000);

end
s.forcePeriodicCoilVoltages;
% plot results
s.plotCoilCurrents;
s.plotCoilFluxLinkages;
s.plotCoilInducedVoltages;

id_srf = (s.coils.phaseA.current-0.5*s.coils.phaseB.current-0.5*s.coils.phaseC.current)*2/3;
iq_srf = (s.coils.phaseB.current-s.coils.phaseC.current)/sqrt(3);

fld_srf = (s.coils.phaseA.fluxLinkage-0.5*s.coils.phaseB.fluxLinkage-0.5*s.coils.phaseC.fluxLinkage)*2/3;
flq_srf = (s.coils.phaseB.fluxLinkage-s.coils.phaseC.fluxLinkage)/sqrt(3);

Te_dq = (3/2) * (gv_p/2) * (fld_srf.*iq_srf - flq_srf.*id_srf);

figure;
hold on; box on;
plot(s.simTime, te);
plot(s.simTime, te_mst);
plot(s.simTime, Te_dq);
xlabel('Time [s]'); ylabel('Electric torque [Nm]');
legend('Arkkio', 'MST', 'Te_{dq}');
set(gca,'ylim',[0,1.2*max(te)])
hold off;



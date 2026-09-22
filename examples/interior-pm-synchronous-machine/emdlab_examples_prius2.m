%{
note: calculation of the cogging torque of Prius 2
%}

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% design variables
gv_ISD = 185.2;
gv_OSD = 282;
gv_Lstk = 75;
gv_Ns = 48;
gv_p = 8;
gv_wst = 6.5;
gv_dss = 30;
gv_bs0 = 2.4;
gv_hs0 = 1;
gv_tta = 40;
gv_Dsh = 112;
gv_dm = 5;
gv_gap = 0.6;
gv_alpha_v = 160;
gv_Hc = -867100;
gv_wtrib = 0.8;
gv_wrrib = 0;
gv_bm0 = 13.6;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_ipm5(g,gv_Dsh,gv_ISD-2*gv_gap,gv_p,gv_dm,gv_alpha_v,gv_wtrib,gv_wrrib,gv_bm0,'rotor','magnet','rap')

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [6,gv_gap/2,6], r, 'linear','extrap');
g.setMeshLengthByRadialFunction(f_mesh);

% mesh generation
m = g.generateMesh('mg0');

% add materials
m.addMaterial('core', emdlab_mlib_es_M270_35A);
m.addMaterial('copper', emdlab_mlib_copper);

% set materials
m.setMaterial('rotor','core');
m.setMaterial('stator','core');
m.setMaterial('sc','copper');

% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns)
m.aux_cmxjcrj('sap',gv_Ns)
m.aux_cmxjcrj('rotor',gv_p)
m.aux_unify('rap');
m.aux_cmxjcrj('rap',gv_p)
m.aux_cmxjcr('sc',gv_Ns)
m.aux_cmxcr('magnet',gv_p)

% generate air gap mesh
m.aux_addCircularAirGapInterface('ag',0,0,gv_ISD/2-gv_gap,0,0,gv_ISD/2,2)

% getting an instance of solver object
s = emdlab_solvers_mt2d_tl3_ihnlwm(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% set magnetization
tmp = pi/2-deg2rad(gv_alpha_v/2);
u1 = [cos(tmp),-sin(tmp)];
u2 = [u1(1),-u1(2)];
for i = 1:2:gv_p
    u11 = emdlab_g2d_rotatePoints(u1,(i-1)*2*pi/gv_p);
    s.setMagnetization(['magnet1',num2str(i)],gv_Hc,u11);
    m.setMeshZoneColor(['magnet1',num2str(i)],0,125,223);
    u22= emdlab_g2d_rotatePoints(u2,(i-1)*2*pi/gv_p);
    s.setMagnetization(['magnet2',num2str(i)],gv_Hc,u22);
    m.setMeshZoneColor(['magnet2',num2str(i)],0,125,223);
end
u1 = -emdlab_g2d_rotatePoints(u1,2*pi/gv_p);
u2 = -emdlab_g2d_rotatePoints(u2,2*pi/gv_p);
for i = 2:2:gv_p
    u11 = emdlab_g2d_rotatePoints(u1,(i-2)*2*pi/gv_p);
    s.setMagnetization(['magnet1',num2str(i)],gv_Hc,u11);
    m.setMeshZoneColor(['magnet1',num2str(i)],255,70,70);
    u22= emdlab_g2d_rotatePoints(u2,(i-2)*2*pi/gv_p);
    s.setMagnetization(['magnet2',num2str(i)],gv_Hc,u22);
    m.setMeshZoneColor(['magnet2',num2str(i)],255,70,70);
end

% define moving region
s.defineMovingRegion('moving_1', ["rotor", "rap", "magnet1"+string(1:gv_p), "magnet2"+string(1:gv_p)], 'ag');

% apply boundary conditions: zero vector potential on all boundary nodes
s.setAzBC(m.getfbn, 0);

% disable solver monitor
s.setMonitor(0);

% synchronous speed
wm = 6000*pi/30;

% total time steps for 20e-3 simulation
Nt = 40; 

% simulation time step size
timeStep = (deg2rad(7.5)/wm)/Nt;
te = zeros(1,Nt);
s.solveForInitialConditions(20);
for i = 1:Nt    
    % calculate torque
    te(i) = s.evalTorqueByArkkio('ag', gv_gap);
    % rotate moving region
    s.rotateMovingRegion('moving_1', wm*timeStep);
    s.solveForOneTimeStep(timeStep);
end

% cogging torque
figure;
hold on; box on;
rp = rad2deg(s.simTime(1:end-1)*wm); % rotor positions
rp_fine = linspace(rp(1),rp(end),100*Nt);
plot(rp_fine, interp1(rp,te,rp_fine,'spline'));
plot(rp,te,'o');
xlabel('Rotor Position [deg]'); ylabel('Cogging Torque [Nm]');

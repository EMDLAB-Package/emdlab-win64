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
gv_Is = 2.3;
gv_embrace = 0.82;
gv_Hc = -922100;
Nagl = 2;

% define geometry data base
g = emdlab_g2d_db;

% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_spm1(g, gv_Dsh, gv_ISD-2*gv_g, gv_p, gv_dm, gv_embrace, 'rotor', 'magnet', 'rap')

% setting the wireframe mesh by mesh size function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_dm-gv_g,gv_ISD/2,gv_OSD/2], [3,gv_g/Nagl,gv_g/Nagl,4], r, 'linear','extrap');
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
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,Nagl)

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
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
s.defineCoil('phaseA');
s.defineCoil('phaseB');
s.defineCoil('phaseC');
% index of each coil assocciating with each winding
pAn = [1,2,3]; pAn = [pAn, pAn+18];
pAp = [10,11,12]; pAp = [pAp, pAp+18];
pBp = pAp + 6;
pBn = pAn + 6;
pCp = pAn + 3;
pCn = pAp + 3;
% assingation of mesh zones to windings
s.addMeshZone2Coil('phaseA', 'sc'+string(pAp), gv_Ntc, 1);
s.addMeshZone2Coil('phaseA', 'sc'+string(pAn), gv_Ntc, -1);
s.addMeshZone2Coil('phaseB', 'sc'+string(pBp), gv_Ntc, 1);
s.addMeshZone2Coil('phaseB', 'sc'+string(pBn), gv_Ntc, -1);
s.addMeshZone2Coil('phaseC', 'sc'+string(pCp), gv_Ntc, 1);
s.addMeshZone2Coil('phaseC', 'sc'+string(pCn), gv_Ntc, -1);

% apply boundary conditions
s.setAzBC(m.getfbn, 0);

% phase currents in stator reference frame
alpha_d4A = -110 * pi/180;
alpha_Is_srf = @(alpha_Is_rrf) alpha_Is_rrf + alpha_d4A;
id_srf = @(alpha_Is_rrf) gv_Is * cos(alpha_Is_srf(alpha_Is_rrf));
iq_srf = @(alpha_Is_rrf) gv_Is * sin(alpha_Is_srf(alpha_Is_rrf));
ia_srf = @(alpha_Is_rrf) id_srf(alpha_Is_rrf);
ib_srf = @(alpha_Is_rrf) -id_srf(alpha_Is_rrf)/2 + sqrt(3)*iq_srf(alpha_Is_rrf)/2;
ic_srf = @(alpha_Is_rrf) -id_srf(alpha_Is_rrf)/2 - sqrt(3)*iq_srf(alpha_Is_rrf)/2;

% total simulation points for changing the current angle
Nt = 20;

% allocate memory to store calculated electric torque
alpha_Is_rrf = linspace(0,pi,Nt);
te = zeros(1,Nt);
te_mst = zeros(1,Nt);
for i = 1:Nt
    % set coil currents
    s.setCoilCurrent('phaseA', ia_srf(alpha_Is_rrf(i)));
    s.setCoilCurrent('phaseB', ib_srf(alpha_Is_rrf(i)));
    s.setCoilCurrent('phaseC', ic_srf(alpha_Is_rrf(i)));
    % run solver
    s.solve;
    % calculate torque
    te(i) = s.evalTorqueByArkkio('ag', gv_g);
    te_mst(i) = s.evalTorqueByMST3(0,0,gv_ISD/2-gv_g/2,gv_g,1000);
end

% plot static torque curve
figure;
hold on; box on;
plot(alpha_Is_rrf*180/pi, te);
plot(alpha_Is_rrf*180/pi, te_mst);
xlabel('Current angle [deg]'); ylabel('Electric torque [Nm]');
legend('Arkkio', 'MST');
set(gca,'xlim',[0,180]);
set(gca,'ylim',[0,1.1*max(te)]);
hold off;


function [hhcl, p3] = emdlab_m3d_extrudeQM2HHM(qcl, p2, zLevels)
% EMDLAB_M3D_EXTRUDEQM2HHM Extrude 2D quad mesh into 3D hex mesh
%
% qcl     : Nq x 4 quad connectivity
% p2      : Np x 2 or Np x 3 node coordinates
% zLevels : vector of z coordinates, length = nLayers+1
%
% p3      : 3D node coordinates
% hhcl    : (Nq*nLayers) x 8 hexahedral connectivity

if size(p2,2) == 2
    p2 = [p2, zeros(size(p2,1),1)];
elseif size(p2,2) ~= 3
    error('p2 must be N x 2 or N x 3');
end

Np = size(p2,1);
Nq = size(qcl,1);
Nez = numel(zLevels) - 1;

if Nez < 1
    error('zLevels must contain at least two z coordinates');
end

% Build extruded node coordinates
p3 = zeros(Np*(Nez+1), 3);
for k = 1:(Nez+1)
    idx = (1:Np) + (k-1)*Np;
    p3(idx,1) = p2(:,1);
    p3(idx,2) = p2(:,2);
    p3(idx,3) = zLevels(k);
end

% Build hexahedral connectivity
hhcl = zeros(Nq*Nez, 8);
row = 1;

for k = 1:Nez
    offsetBot = (k-1)*Np;
    offsetTop = k*Np;

    bot = qcl + offsetBot;
    top = qcl + offsetTop;

    % Reverse top face ordering for consistent hex orientation
    hhcl(row:row+Nq-1,:) = [bot(:,[1,2,3,4]), top(:,[1,2,3,4])];

    row = row + Nq;
end
end

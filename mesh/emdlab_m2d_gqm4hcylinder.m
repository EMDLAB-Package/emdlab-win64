function meshZone = emdlab_m2d_gqm4hcylinder(x0,y0,r1,r2,meshSize)
%EMDLAB_M2D_GQM4HOLLOWCYLINDER Generate a structured quad mesh for an annulus
%
%   meshZone = emdlab_m2d_gqm4hollowCylinder(x0,y0,r1,r2,meshSize)
%
% Inputs
%   x0, y0   : center
%   r1       : inner radius
%   r2       : outer radius
%   meshSize : target element size
%
% Output
%   meshZone : mesh object created by emdlab_m2d_qmz(elem,node)

if meshSize <= 0
    error('meshSize must be positive.');
end
if r1 <= 0 || r2 <= 0
    error('r1 and r2 must be positive.');
end
if r2 <= r1
    error('Require r2 > r1.');
end

% Circumferential divisions from the mean radius.
rMid = 0.5 * (r1 + r2);
nTheta = max(8, ceil(2*pi*r1 / meshSize));
nTheta = 4 * ceil(nTheta / 4);

% Radial divisions from thickness.
nRad = max(1, ceil((r2 - r1) / meshSize));

theta = linspace(0, 2*pi, nTheta + 1);
rho = linspace(r1, r2, nRad + 1);

[T, R] = meshgrid(theta, rho);

X = x0 + R .* cos(T);
Y = y0 + R .* sin(T);

node = zeros((nRad + 1) * nTheta, 2);
id = zeros(nRad + 1, nTheta);

k = 0;
for j = 1:(nRad + 1)
    for i = 1:nTheta
        k = k + 1;
        node(k,:) = [X(j,i), Y(j,i)];
        id(j,i) = k;
    end
end

elem = zeros(nRad * nTheta, 4);
e = 0;
for j = 1:nRad
    for i = 1:nTheta
        i2 = i + 1;
        if i == nTheta
            i2 = 1;
        end

        e = e + 1;
        elem(e,:) = [id(j,i), id(j,i2), id(j+1,i2), id(j+1,i)];
    end
end

meshZone = emdlab_m2d_qmz(elem, node);
end

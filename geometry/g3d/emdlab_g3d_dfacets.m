% inputs:
% 1: facets
% 2: vertices
% 3: query points

% output is signed distances

function d = emdlab_g3d_dfacets(f, v, p)

% number of points
Np = size(p,1);

% x and y coordinates of the start and end point of the segment
% segments are from point a to b
xa = v(f(:,1),1);
ya = v(f(:,1),2);
xb = v(f(:,2),1);
yb = v(f(:,2),2);

u = [xb - xa, yb - ya];
uMAG = vecnorm(u,2,2);
u = u./uMAG;

% signed distances output
d = zeros(Np,1);

for i = 1:Np

    xp = p(i,1);
    yp = p(i,2);

    ap = [xp-xa,yp-ya];
    apMAG = vecnorm(ap,2,2);

    bp = [xp-xb,yp-yb];
    bpMAG = vecnorm(bp,2,2);

    % dot product
    alpha = sum(ap.*u,2);

    % cross product
    beta = sign(sum(ap(:,1).*u(:,2)-ap(:,2).*u(:,1),2));

    y = vecnorm(ap-alpha.*u,2,2);

    index = alpha<0;
    y(index) = apMAG(index);

    index = alpha>uMAG;
    y(index) = bpMAG(index);

    [d(i),index] = min(y);
    d(i) = d(i) * beta(index);

end

end
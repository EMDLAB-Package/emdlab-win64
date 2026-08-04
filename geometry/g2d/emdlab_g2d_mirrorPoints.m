%{
EMDLAB: Electrical Machines Design Laboratory

=> mirror points on x-y plane
=> p is a matrix [Np x 2]

This function mirrors points with respect to the p0-p1 line.

=> p0 = [x0, y0]
=> p1 = [x1, y1]

Note:
When x1 and y1 are not defined, they are considered equal to 0.
Therefore, the mirror line passes through the origin and [x0, y0].
%}

function newP = emdlab_g2d_mirrorPoints(p, x0, y0, x1, y1)

if nargin < 4
    x1 = 0;
    y1 = 0;
end

if size(p, 2) ~= 2
    error('The size of point matrix must be [Npx2].');
end

% Direction vector of the mirror line
u = [x1 - x0, y1 - y0];

% Check that p0 and p1 are not identical
nu = norm(u);

if nu == 0
    error('The two points defining the mirror line must be different.');
end

% Unit direction vector
u = u / nu;

% Move points so that p0 becomes the origin
q = p - [x0, y0];

% Projection of each point onto the mirror line
qproj = sum(q .* u, 2) .* u;

% Mirror the points
% q_mirror = 2*qproj - q
qmirror = 2 * qproj - q;

% Move points back to the original coordinate system
newP = qmirror + [x0, y0];

end

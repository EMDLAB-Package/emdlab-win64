%{ 
EMDLAB: Electrical Machines Design Laboratory

=> mirror points on x-y plane
=> p is a matrix [Np x 2]

this function mirror points with respect to p0p1 line
=> p0 = [x0,y0]
=> p1 = [x1,y1]

Note: when x1 and y1 are not defined they will be considered equal to 0 so
the line is passing through origin and x0 and y0

%}

function [newX, newY] = emdlab_g2d_mirrorPointsXY(x, y, x0, y0, x1, y1)

if nargin < 5
    x1 = 0;
    y1 = 0;
end

if ~all(size(x) == size(y))
    error('x and y vectors must have the same size.');
end

% Direction vector of the mirror line
u = [x1 - x0, y1 - y0];

if u(1) == 0 && u(2) == 0
    error('The mirror line is not valid because the two line points are identical.');
end

% Normalize direction vector
u = u/norm(u);

% Shift points so line passes through origin
x = x - x0;
y = y - y0;

% Projection onto line
proj = (x*ux + y*uy).*u;

% Perpendicular component
perp = p - proj;

% Shift back
newP = proj - perp + [x0,y0];

end
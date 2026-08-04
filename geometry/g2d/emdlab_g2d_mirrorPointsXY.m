%{
EMDLAB: Electrical Machines Design Laboratory

=> mirror points on x-y plane
=> x and y contain the point coordinates

This function mirrors points with respect to the p0-p1 line.

=> p0 = [x0, y0]
=> p1 = [x1, y1]

Note:
When x1 and y1 are not defined, they are considered equal to 0.
Therefore, the line passes through the origin and [x0, y0].
%}

function [newX, newY] = emdlab_g2d_mirrorPointsXY(x, y, x0, y0, x1, y1)

if nargin < 5
    x1 = 0;
    y1 = 0;
end

if ~isequal(size(x), size(y))
    error('x and y vectors must have the same size.');
end

% Direction vector of the mirror line
ux = x1 - x0;
uy = y1 - y0;

% Check validity of the mirror line
normU = sqrt(ux^2 + uy^2);

if normU == 0
    error('The mirror line is not valid because the two line points are identical.');
end

% Normalize direction vector
ux = ux / normU;
uy = uy / normU;

% Shift points so that p0 becomes the origin
x = x - x0;
y = y - y0;

% Projection of each point onto the mirror line
proj = x .* ux + y .* uy;

% Mirrored coordinates
newX = 2 .* proj .* ux - x;
newY = 2 .* proj .* uy - y;

% Shift points back to the original coordinate system
newX = newX + x0;
newY = newY + y0;

end

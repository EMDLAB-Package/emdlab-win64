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

function newP = emdlab_g2d_mirrorPoints(p, x0, y0, x1, y1)

if nargin < 4
    x1 = 0;
    y1 = 0;
end

if size(p,2) ~= 2
    error('The size of point matrix must be [Npx2].');
end

u = [x1 - x0, y1 - y0];
u = u/norm(u);

p = p - [x0,y0];
p - sum(p.*u,2).*u
newP = zeros(size(p,1),2);

p(:,1) = p(:,1) - xc;
p(:,2) = p(:,2) - yc;

newP(:,1) = p(:,1) * cos(rotAngle) - p(:,2) * sin(rotAngle);
newP(:,2) = p(:,1) * sin(rotAngle) + p(:,2) * cos(rotAngle);

newP(:,1) = newP(:,1) + xc;
newP(:,2) = newP(:,2) + yc;

end
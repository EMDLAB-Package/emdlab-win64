% EMDLAB: Electrical Machines Design Laboratory
% shift points on x-y plane
% p is a matrix [Np x 2]
% dx, dy: amount of shift in x- and y-directions

function newP = emdlab_g2d_shiftPoints(p, dx, dy)

if size(p,2) ~= 2
    error('The size of point matrix must be [Npx2].');
end

newP = p;

newP(:,1) = newP(:,1) + dx;
newP(:,2) = newP(:,2) + dy;

end
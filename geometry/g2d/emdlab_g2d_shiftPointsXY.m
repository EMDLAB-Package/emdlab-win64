% EMDLAB: Electrical Machines Design Laboratory
% shift points on x-y plane
% dx, dy: amount of shift in x- and y-directions

function [xs, ys] = emdlab_g2d_shiftPointsXY(x, y, dx, dy)

xs = x + dx;
ys = y + dy;

end
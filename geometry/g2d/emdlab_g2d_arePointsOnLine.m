function [tf, sortedPts] = emdlab_g2d_arePointsOnLine(pts, x0, y0, ux, uy, tol)

    % pts      : N x 2 points
    % x0, y0   : reference point on the line
    % ux, uy   : line direction vector
    % tol      : distance tolerance
    %
    % tf       : true if all points are on the line
    % sortedPts: points sorted along the direction [ux, uy]

    if nargin < 6
        tol = 1e-6;
    end

    % Direction vector magnitude
    normU = hypot(ux, uy);

    if normU == 0
        error('The line direction vector [ux, uy] cannot be zero.');
    end

    % Vector from (x0,y0) to each point
    dx = pts(:,1) - x0;
    dy = pts(:,2) - y0;

    % 2D cross product
    crossVal = dx .* uy - dy .* ux;

    % Perpendicular distance from each point to the line
    distance = abs(crossVal) ./ normU;

    % Check whether all points lie on the line
    tf = all(distance <= tol);

    % Scalar projection along direction vector [ux, uy]
    projection = dx .* ux + dy .* uy;

    % Sort along the direction vector
    [~, index] = sort(projection);

    sortedPts = pts(index,:);

end
function [points, index] = emdlab_g2d_sortPointsCCW(points, center)

if nargin == 1
    center = [0,0];
end

% angles
a = atan_02pi([points(:, 1) - center(1), points(:, 2) - center(2)]);

[~, index] = sort(a);
points = points(index,:);

while true
    
    if points(1,1)*points(end,2) -  points(end,1)*points(1,2) > -1e-6
        break;
    end
    
    points = circshift(points,1);
    index = circshift(index,1);
    
end

end
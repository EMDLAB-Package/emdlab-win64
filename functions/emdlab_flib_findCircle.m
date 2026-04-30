function [xc,yc,r] = emdlab_flib_findCircle(x1,y1,x2,y2,x3,y3)

ab = [x2-x1,y2-y1]; 
bc = [x3-x2,y3-y2];

ab = ab/norm(ab);
bc = bc/norm(bc);

[ux1,uy1] = rotate_pts(ab(1),ab(2),pi/2);
[ux2,uy2] = rotate_pts(bc(1),bc(2),pi/2);

[xc,yc] = getLineLineIntersection([x1+x2,y1+y2]/2, [ux1,uy1], [x2+x3,y2+y3]/2, [ux2,uy2]);

r = norm([x1-xc, y1-yc]);


end
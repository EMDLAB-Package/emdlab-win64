function [m1, m2] = emdlab_flib_addSlidingCircularAirGapInterface(m, ipts, opts, tol)

% sort points

r1 = mean(vecnorm(ipts,2,2));
r2 = mean(vecnorm(opts,2,2));

rmid = (r1+r2)/2;



end
function [B_ext,nu_ext] = emdlab_extend_linear_to_air(B,nu)

mu0 = 4*pi*1e-7;
nu0 = 1/mu0;

% derivative at last point
pp = pchip(B,nu);
breaks = pp.breaks;
coefs  = pp.coefs;   % Each row is [a b c d]
% Differentiate polynomial coefficients
% Original: a*(x-xk)^3 + b*(x-xk)^2 + c*(x-xk) + d
% Derivative: 3a*(x-xk)^2 + 2b*(x-xk) + c
coefs_der = [3*coefs(:,1), 2*coefs(:,2), coefs(:,3)];
pp_der = mkpp(breaks, coefs_der);
slope = ppval(pp_der, B(end));

% distance needed to reach air reluctivity
DeltaB = (nu0 - nu(end)) / slope;

% extend moderately
B_ext = [B;
         B(end) + DeltaB;
         B(end) + 1000*DeltaB];

nu_ext = [nu;
          nu0;
          nu0];

end

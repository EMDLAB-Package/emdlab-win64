function pgw = emdlab_flib_propyleneGlycolWater6040Properties(Tavg)
%EMDLAB_FLIB_PROPYLENEGLYCOLWATER6040PROPERTIES
% Thermophysical properties of 60/40 volume% propylene glycol-water mixture.
%
% Input
%   Tavg : Average coolant temperature [degC]
%
% Output
%   pgw.T_C
%   pgw.T_K
%   pgw.p
%   pgw.rho
%   pgw.Cp
%   pgw.k
%   pgw.mu
%   pgw.nu
%   pgw.alpha
%   pgw.Pr
%   pgw.beta
%
% Valid range:
%   0 <= Tavg <= 120 degC
%
% Properties correspond to approximately:
%   60% propylene glycol + 40% water by volume.
%
% The mixture is assumed to remain in liquid state.

%% Input checking

validateattributes(Tavg,{'numeric'},...
    {'real','finite','nonempty'},mfilename,'Tavg');

if any(Tavg(:)<0 | Tavg(:)>120)
    error('Tavg must be between 0 and 120 degC.');
end

%% Temperature

pgw.T_C = Tavg;
pgw.T_K = Tavg + 273.15;

pgw.p = 101325;

%% Temperature table

Ttab = (0:10:120)';

%% Density [kg/m^3]

rho_tab = [
1068;
1063;
1058;
1053;
1048;
1043;
1038;
1033;
1028;
1023;
1018;
1013;
1008];

%% Specific heat [J/(kg.K)]

Cp_tab = [
3100;
3150;
3200;
3250;
3300;
3350;
3400;
3450;
3500;
3550;
3600;
3650;
3700];

%% Thermal conductivity [W/(m.K)]

k_tab = [
0.355;
0.360;
0.365;
0.370;
0.375;
0.380;
0.385;
0.390;
0.395;
0.400;
0.405;
0.410;
0.415];

%% Dynamic viscosity [Pa.s]

mu_tab = [
0.0120;
0.0080;
0.0055;
0.0039;
0.0029;
0.00225;
0.00175;
0.00140;
0.00115;
0.00095;
0.00080;
0.00068;
0.00058];

%% Prandtl number [-]

Pr_tab = [
105;
72;
51;
38;
29;
23;
18;
14.5;
12;
10;
8.5;
7.2;
6.2];

%% Thermal expansion coefficient [1/K]

beta_tab = [
5.2e-4;
5.5e-4;
5.8e-4;
6.1e-4;
6.4e-4;
6.7e-4;
7.0e-4;
7.3e-4;
7.6e-4;
7.9e-4;
8.2e-4;
8.5e-4;
8.8e-4];

%% Interpolation

pgw.rho = interp1(Ttab,rho_tab,Tavg,'pchip');

pgw.Cp = interp1(Ttab,Cp_tab,Tavg,'pchip');

pgw.k = interp1(Ttab,k_tab,Tavg,'pchip');

pgw.mu = interp1(Ttab,mu_tab,Tavg,'pchip');

pgw.Pr = interp1(Ttab,Pr_tab,Tavg,'pchip');

pgw.beta = interp1(Ttab,beta_tab,Tavg,'pchip');

%% Derived properties

pgw.nu = pgw.mu ./ pgw.rho;

pgw.alpha = pgw.k ./ (pgw.rho .* pgw.Cp);

end
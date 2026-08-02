function oil = emdlab_flib_motorJetOilProperties(Tavg)
%EMDLAB_FLIB_MOTORJETOILPROPERTIES Thermophysical properties of
%synthetic motor jet oil.
%
% Input
%   Tavg : Average oil temperature [degC]
%
% Output
%   oil.T_C
%   oil.T_K
%   oil.p
%   oil.rho
%   oil.Cp
%   oil.k
%   oil.mu
%   oil.nu
%   oil.alpha
%   oil.Pr
%   oil.beta
%
% Valid range:
%   0 <= Tavg <= 180 degC
%
% Properties represent a typical synthetic jet oil used for
% high-temperature motor cooling applications.
%
% The oil is assumed to remain in liquid state.

%% Input checking

validateattributes(Tavg,{'numeric'},...
    {'real','finite','nonempty'},mfilename,'Tavg');

if any(Tavg(:)<0 | Tavg(:)>180)
    error('Tavg must be between 0 and 180 degC.');
end

%% Temperature

oil.T_C = Tavg;
oil.T_K = Tavg + 273.15;

oil.p = 101325;

%% Temperature table

Ttab = (0:10:180)';

%% Density [kg/m^3]

rho_tab = [
980;
976;
972;
968;
964;
960;
956;
952;
948;
944;
940;
936;
932;
928;
924;
920;
916;
912;
908];

%% Specific heat [J/(kg.K)]

Cp_tab = [
1950;
1980;
2010;
2040;
2070;
2100;
2130;
2160;
2190;
2220;
2250;
2280;
2310;
2340;
2370;
2400;
2430;
2460;
2490];

%% Thermal conductivity [W/(m.K)]

k_tab = [
0.140;
0.139;
0.138;
0.137;
0.136;
0.135;
0.134;
0.133;
0.132;
0.131;
0.130;
0.129;
0.128;
0.127;
0.126;
0.125;
0.124;
0.123;
0.122];

%% Dynamic viscosity [Pa.s]

mu_tab = [
0.032;
0.022;
0.0155;
0.0115;
0.0088;
0.0069;
0.0055;
0.0045;
0.0037;
0.0031;
0.0027;
0.0023;
0.0020;
0.00175;
0.00155;
0.00138;
0.00125;
0.00113;
0.00103];

%% Prandtl number [-]

Pr_tab = [
445;
314;
226;
171;
134;
107;
88;
73;
61;
53;
47;
41;
36;
32;
29;
26;
24;
22;
21];

%% Thermal expansion coefficient [1/K]

beta_tab = [
7.0e-4;
7.2e-4;
7.4e-4;
7.6e-4;
7.8e-4;
8.0e-4;
8.2e-4;
8.4e-4;
8.6e-4;
8.8e-4;
9.0e-4;
9.2e-4;
9.4e-4;
9.6e-4;
9.8e-4;
1.0e-3;
1.02e-3;
1.04e-3;
1.06e-3];

%% Interpolation

oil.rho = interp1(Ttab,rho_tab,Tavg,'pchip');

oil.Cp = interp1(Ttab,Cp_tab,Tavg,'pchip');

oil.k = interp1(Ttab,k_tab,Tavg,'pchip');

oil.mu = interp1(Ttab,mu_tab,Tavg,'pchip');

oil.Pr = interp1(Ttab,Pr_tab,Tavg,'pchip');

oil.beta = interp1(Ttab,beta_tab,Tavg,'pchip');

%% Derived properties

oil.nu = oil.mu ./ oil.rho;

oil.alpha = oil.k ./ (oil.rho .* oil.Cp);

end
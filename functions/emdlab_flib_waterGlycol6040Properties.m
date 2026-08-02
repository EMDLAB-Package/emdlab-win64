function wg = emdlab_flib_waterGlycol6040Properties(Tavg)
%EMDLAB_FLIB_WATERGLYCOL6040PROPERTIES Thermophysical properties of
%60/40 volume% water-ethylene glycol mixture.
%
% Input
%   Tavg : Average coolant temperature [degC]
%
% Output
%   wg.T_C
%   wg.T_K
%   wg.p
%   wg.rho
%   wg.Cp
%   wg.k
%   wg.mu
%   wg.nu
%   wg.alpha
%   wg.Pr
%   wg.beta
%
% Valid range:
%   0 <= Tavg <= 120 degC
%
% Properties correspond to approximately 60/40 volume %
% water-ethylene glycol mixture.
%
% The mixture is assumed to remain in liquid state.

%% Input checking

validateattributes(Tavg,{'numeric'},...
    {'real','finite','nonempty'},mfilename,'Tavg');

if any(Tavg(:)<0 | Tavg(:)>120)
    error('Tavg must be between 0 and 120 degC.');
end

%% Temperature

wg.T_C = Tavg;
wg.T_K = Tavg + 273.15;

wg.p = 101325;

%% Temperature table

Ttab = (0:10:120)';

%% Density [kg/m^3]

rho_tab = [
1065;
1060;
1054;
1048;
1042;
1036;
1030;
1024;
1018;
1012;
1006;
1000;
994];

%% Specific heat [J/(kg.K)]

Cp_tab = [
3500;
3560;
3620;
3680;
3740;
3800;
3860;
3920;
3980;
4040;
4100;
4160;
4220];

%% Thermal conductivity [W/(m.K)]

k_tab = [
0.405;
0.410;
0.415;
0.420;
0.425;
0.430;
0.435;
0.440;
0.445;
0.450;
0.455;
0.460;
0.465];

%% Dynamic viscosity [Pa.s]

mu_tab = [
3.80e-3;
2.60e-3;
1.85e-3;
1.38e-3;
1.05e-3;
8.30e-4;
6.70e-4;
5.60e-4;
4.80e-4;
4.20e-4;
3.70e-4;
3.30e-4;
3.00e-4];

%% Prandtl number [-]

Pr_tab = [
33;
23;
17;
13;
10.5;
8.8;
7.5;
6.5;
5.7;
5.1;
4.6;
4.2;
3.9];

%% Thermal expansion coefficient [1/K]

beta_tab = [
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
8.8e-4;
9.1e-4];

%% Interpolation

wg.rho  = interp1(Ttab,rho_tab,Tavg,'pchip');

wg.Cp   = interp1(Ttab,Cp_tab,Tavg,'pchip');

wg.k    = interp1(Ttab,k_tab,Tavg,'pchip');

wg.mu   = interp1(Ttab,mu_tab,Tavg,'pchip');

wg.Pr   = interp1(Ttab,Pr_tab,Tavg,'pchip');

wg.beta = interp1(Ttab,beta_tab,Tavg,'pchip');

%% Derived properties

wg.nu = wg.mu ./ wg.rho;

wg.alpha = wg.k ./ (wg.rho .* wg.Cp);

end
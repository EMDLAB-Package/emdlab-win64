function wg = emdlab_flib_waterGlycol5050Properties(Tavg)
%EMDLAB_FLIB_WATERGLYCOL5050PROPERTIES Thermophysical properties of
%50/50 volume% water-ethylene glycol mixture.
%
% Input
%   Tavg : Average coolant temperature [degC]
%
% Output
%   wg.T_C     Temperature [degC]
%   wg.T_K     Temperature [K]
%   wg.p       Reference pressure [Pa]
%   wg.rho     Density [kg/m^3]
%   wg.Cp      Specific heat capacity [J/(kg.K)]
%   wg.k       Thermal conductivity [W/(m.K)]
%   wg.mu      Dynamic viscosity [Pa.s]
%   wg.nu      Kinematic viscosity [m^2/s]
%   wg.alpha   Thermal diffusivity [m^2/s]
%   wg.Pr      Prandtl number [-]
%   wg.beta    Volumetric thermal expansion coefficient [1/K]
%
% Valid range:
%   0 <= Tavg <= 120 degC
%
% Properties correspond to approximately 50/50 volume %
% ethylene glycol-water coolant.
%
% Note:
%   The mixture is assumed to remain in liquid state.
%   For automotive/electrical cooling systems, the actual pressure
%   is usually higher than atmospheric pressure.

%% Input checking

validateattributes(Tavg,{'numeric'},...
    {'real','finite','nonempty'},mfilename,'Tavg');

if any(Tavg(:)<0 | Tavg(:)>120)
    error('Tavg must be between 0 and 120 degC.');
end

%% Temperature

wg.T_C = Tavg;
wg.T_K = Tavg + 273.15;

% Reference pressure
wg.p = 101325;

%% Temperature table

Ttab = (0:10:120)';

%% Density [kg/m^3]

rho_tab = [
1088;
1083;
1078;
1072;
1066;
1060;
1054;
1048;
1041;
1035;
1028;
1021;
1014];

%% Specific heat [J/(kg.K)]

Cp_tab = [
3300;
3370;
3450;
3530;
3610;
3690;
3770;
3850;
3930;
4010;
4090;
4170;
4250];

%% Thermal conductivity [W/(m.K)]

k_tab = [
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
0.415;
0.420];

%% Dynamic viscosity [Pa.s]

mu_tab = [
6.00e-3;
4.00e-3;
2.90e-3;
2.20e-3;
1.75e-3;
1.40e-3;
1.15e-3;
9.70e-4;
8.30e-4;
7.20e-4;
6.30e-4;
5.60e-4;
5.00e-4];

%% Prandtl number [-]

Pr_tab = [
55;
37;
27;
21;
17;
13.5;
11;
9.0;
7.5;
6.5;
5.7;
5.1;
4.6];

%% Thermal expansion coefficient [1/K]

beta_tab = [
5.0e-4;
5.3e-4;
5.6e-4;
5.9e-4;
6.2e-4;
6.5e-4;
6.8e-4;
7.1e-4;
7.4e-4;
7.7e-4;
8.0e-4;
8.3e-4;
8.6e-4];

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
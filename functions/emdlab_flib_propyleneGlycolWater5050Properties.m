function pgw = emdlab_flib_propyleneGlycolWater5050Properties(Tavg)
%EMDLAB_FLIB_PROPYLENEGLYCOLWATER5050PROPERTIES
% Thermophysical properties of 50/50 volume% propylene glycol-water mixture.
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
%   50% propylene glycol + 50% water by volume.
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
1065;
1060;
1055;
1050;
1045;
1040;
1035;
1030;
1025;
1020;
1015;
1010;
1005];

%% Specific heat [J/(kg.K)]

Cp_tab = [
3400;
3440;
3490;
3540;
3590;
3640;
3690;
3740;
3790;
3840;
3890;
3940;
3990];

%% Thermal conductivity [W/(m.K)]

k_tab = [
0.380;
0.385;
0.390;
0.395;
0.400;
0.405;
0.410;
0.415;
0.420;
0.425;
0.430;
0.435;
0.440];

%% Dynamic viscosity [Pa.s]

mu_tab = [
0.0090;
0.0060;
0.0042;
0.0030;
0.0022;
0.0017;
0.00135;
0.00110;
0.00090;
0.00076;
0.00065;
0.00057;
0.00050];

%% Prandtl number [-]

Pr_tab = [
80;
54;
38;
28;
22;
18;
15;
12.5;
10.5;
9.0;
7.8;
6.8;
6.0];

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
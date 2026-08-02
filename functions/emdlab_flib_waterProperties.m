function water = emdlab_flib_waterProperties(Tavg)
%WATERPROPERTIES Thermophysical properties of liquid water.
%
% Input
%   Tavg : Average water temperature [°C]
%
% Output
%   water.T_C
%   water.T_K
%   water.p
%   water.rho
%   water.Cp
%   water.k
%   water.mu
%   water.nu
%   water.alpha
%   water.Pr
%   water.beta
%
% Valid for:
%   0 <= Tavg <= 200 °C
%
% Properties correspond to liquid water (sufficient pressure to
% remain liquid).

%% Input checking

validateattributes(Tavg,{'numeric'},...
    {'real','finite','nonempty'},mfilename,'Tavg');

if any(Tavg(:)<0 | Tavg(:)>200)
    error('Tavg must be between 0 and 200 °C.');
end

%% Temperature

water.T_C = Tavg;
water.T_K = Tavg + 273.15;

% Nominal pressure
water.p = 101325;

%% Temperature table

Ttab = (0:10:200)';

%% Density [kg/m^3]

rho_tab = [...
999.84;
999.70;
998.21;
995.65;
992.22;
988.04;
983.20;
977.76;
971.80;
965.31;
958.40;
951.00;
943.20;
934.90;
926.10;
916.80;
907.00;
896.80;
886.10;
875.00;
863.40];

%% Specific heat [J/(kg.K)]

Cp_tab = [...
4217;
4195;
4182;
4179;
4179;
4181;
4184;
4188;
4193;
4198;
4205;
4212;
4219;
4227;
4235;
4243;
4252;
4261;
4270;
4280;
4290];

%% Thermal conductivity [W/(m.K)]

k_tab = [...
0.561;
0.579;
0.598;
0.615;
0.630;
0.643;
0.653;
0.661;
0.667;
0.671;
0.674;
0.675;
0.674;
0.672;
0.668;
0.663;
0.657;
0.650;
0.642;
0.633;
0.624];

%% Dynamic viscosity [Pa.s]

mu_tab = [...
1.792e-3;
1.307e-3;
1.002e-3;
8.01e-4;
6.53e-4;
5.47e-4;
4.66e-4;
4.04e-4;
3.55e-4;
3.16e-4;
2.82e-4;
2.57e-4;
2.37e-4;
2.20e-4;
2.05e-4;
1.92e-4;
1.81e-4;
1.71e-4;
1.62e-4;
1.54e-4;
1.47e-4];

%% Prandtl number

Pr_tab = [...
13.46;
9.45;
7.01;
5.42;
4.31;
3.56;
3.02;
2.62;
2.32;
2.08;
1.90;
1.74;
1.61;
1.50;
1.40;
1.32;
1.25;
1.19;
1.14;
1.09;
1.05];

%% Volumetric expansion coefficient [1/K]

beta_tab = [...
5.0e-5;
8.8e-5;
2.1e-4;
3.0e-4;
3.9e-4;
4.6e-4;
5.2e-4;
5.8e-4;
6.3e-4;
6.8e-4;
7.3e-4;
7.7e-4;
8.1e-4;
8.5e-4;
8.9e-4;
9.3e-4;
9.7e-4;
1.01e-3;
1.05e-3;
1.08e-3;
1.12e-3];

%% Interpolation

water.rho = interp1(Ttab,rho_tab,Tavg,'pchip');
water.Cp  = interp1(Ttab,Cp_tab,Tavg,'pchip');
water.k   = interp1(Ttab,k_tab,Tavg,'pchip');
water.mu  = interp1(Ttab,mu_tab,Tavg,'pchip');
water.Pr  = interp1(Ttab,Pr_tab,Tavg,'pchip');
water.beta = interp1(Ttab,beta_tab,Tavg,'pchip');

%% Derived properties

water.nu = water.mu ./ water.rho;

water.alpha = water.k ./ (water.rho .* water.Cp);

end
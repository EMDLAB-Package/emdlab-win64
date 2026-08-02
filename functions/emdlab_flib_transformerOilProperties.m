function oil = emdlab_flib_transformerOilProperties(Tavg)
%TRANSFORMEROILPROPERTIES Thermophysical properties of mineral transformer oil.
%
% Input
%   Tavg : Average oil temperature [°C]
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
% Valid for:
%   0 <= Tavg <= 200 °C
%
% Properties represent a typical uninhibited mineral transformer oil
% (IEC 60296). Values are representative of commonly used transformer
% oils and are suitable for engineering thermal calculations.

%% Input checking

validateattributes(Tavg,{'numeric'},...
    {'real','finite','nonempty'},mfilename,'Tavg');

if any(Tavg(:)<0 | Tavg(:)>200)
    error('Tavg must be between 0 and 200 °C.');
end

%% Temperature

oil.T_C = Tavg;
oil.T_K = Tavg + 273.15;

oil.p = 101325;

%% Temperature table (°C)

Ttab = (0:10:200)';

%% Density [kg/m^3]

rho_tab = [...
895;
888;
881;
874;
867;
860;
853;
846;
839;
832;
825;
818;
811;
804;
797;
790;
783;
776;
769;
762;
755];

%% Specific heat [J/(kg.K)]

Cp_tab = [...
1800;
1830;
1860;
1890;
1920;
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
2400];

%% Thermal conductivity [W/(m.K)]

k_tab = [...
0.130;
0.129;
0.128;
0.127;
0.126;
0.125;
0.124;
0.123;
0.122;
0.121;
0.120;
0.119;
0.118;
0.117;
0.116;
0.115;
0.114;
0.113;
0.112;
0.111;
0.110];

%% Dynamic viscosity [Pa.s]

mu_tab = [...
0.0850;
0.0500;
0.0310;
0.0205;
0.0145;
0.0105;
0.0080;
0.0062;
0.0050;
0.0041;
0.0034;
0.0029;
0.0025;
0.0022;
0.00195;
0.00175;
0.00158;
0.00143;
0.00131;
0.00120;
0.00111];

%% Prandtl number [-]

Pr_tab = [...
1177;
710;
451;
305;
221;
164;
128;
101;
84;
70;
60;
52;
46;
41;
37;
34;
31;
29;
27;
25;
24];

%% Volumetric expansion coefficient [1/K]

beta_tab = [...
7.0e-4;
7.1e-4;
7.2e-4;
7.3e-4;
7.4e-4;
7.5e-4;
7.6e-4;
7.7e-4;
7.8e-4;
7.9e-4;
8.0e-4;
8.1e-4;
8.2e-4;
8.3e-4;
8.4e-4;
8.5e-4;
8.6e-4;
8.7e-4;
8.8e-4;
8.9e-4;
9.0e-4];

%% Interpolation

oil.rho  = interp1(Ttab,rho_tab,Tavg,'pchip');
oil.Cp   = interp1(Ttab,Cp_tab,Tavg,'pchip');
oil.k    = interp1(Ttab,k_tab,Tavg,'pchip');
oil.mu   = interp1(Ttab,mu_tab,Tavg,'pchip');
oil.Pr   = interp1(Ttab,Pr_tab,Tavg,'pchip');
oil.beta = interp1(Ttab,beta_tab,Tavg,'pchip');

%% Derived properties

oil.nu = oil.mu ./ oil.rho;

oil.alpha = oil.k ./ (oil.rho .* oil.Cp);

end
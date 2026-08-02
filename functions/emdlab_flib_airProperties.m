function air = emdlab_flib_airProperties(Tavg)
%AIRPROPERTIES Thermophysical properties of dry air at 1 atm.
%
% Input
%   Tavg : Average air temperature [°C]
%
% Output
%   air.T_C     Temperature [°C]
%   air.T_K     Temperature [K]
%   air.p       Pressure [Pa]
%   air.rho     Density [kg/m^3]
%   air.Cp      Specific heat [J/(kg.K)]
%   air.k       Thermal conductivity [W/(m.K)]
%   air.mu      Dynamic viscosity [Pa.s]
%   air.nu      Kinematic viscosity [m^2/s]
%   air.alpha   Thermal diffusivity [m^2/s]
%   air.Pr      Prandtl number [-]
%   air.beta    Thermal expansion coefficient [1/K]
%
% Temperature data:
%   -100 to 500 °C with 10 °C increment
%
% Outside range:
%   Properties are extrapolated.

%% Input checking

validateattributes(Tavg,{'numeric'},...
    {'real','finite','nonempty'},mfilename,'Tavg');


%% Temperature

air.T_C = Tavg;
air.T_K = Tavg + 273.15;
air.p   = 101325;


%% Temperature table

Ttab = (-100:10:500)';


%% Density [kg/m^3]
% Ideal gas equation

rho_tab = 101325 ./ (287.05*(Ttab+273.15));


%% Specific heat [J/(kg.K)]
% Smooth approximation for dry air

Cp_tab = 1006 + 0.10*Ttab + 1.5e-4*Ttab.^2;


%% Thermal conductivity [W/(m.K)]
% Polynomial fit

k_tab = 0.0241 + ...
        7.5e-5*Ttab - ...
        1.0e-8*Ttab.^2;


%% Dynamic viscosity [Pa.s]
% Sutherland equation

mu_tab = 1.716e-5 .* ...
    ((Ttab+273.15)/273.15).^1.5 .* ...
    (273.15+111)./(Ttab+273.15+111);


%% Prandtl number [-]

Pr_tab = 0.735 - 8.0e-5*Ttab;


%% Check table consistency

assert(length(Ttab)==length(rho_tab),'rho table length error');
assert(length(Ttab)==length(Cp_tab),'Cp table length error');
assert(length(Ttab)==length(k_tab),'k table length error');
assert(length(Ttab)==length(mu_tab),'mu table length error');
assert(length(Ttab)==length(Pr_tab),'Pr table length error');


%% Interpolation with extrapolation

air.rho = interp1(Ttab,rho_tab,Tavg,'pchip','extrap');

air.Cp  = interp1(Ttab,Cp_tab,Tavg,'pchip','extrap');

air.k   = interp1(Ttab,k_tab,Tavg,'pchip','extrap');

air.mu  = interp1(Ttab,mu_tab,Tavg,'pchip','extrap');

air.Pr  = interp1(Ttab,Pr_tab,Tavg,'pchip','extrap');


%% Derived properties

air.nu = air.mu ./ air.rho;

air.alpha = air.k ./ (air.rho .* air.Cp);

air.beta = 1 ./ air.T_K;


end
function hydrogen = emdlab_flib_hydrogenProperties(Tavg)
%HYDROGENPROPERTIES Thermophysical properties of dry hydrogen at 1 atm.
%
% Input
%   Tavg : Average hydrogen temperature [°C]
%
% Output
%   h2.T_C     Temperature [°C]
%   h2.T_K     Temperature [K]
%   h2.p       Pressure [Pa]
%   h2.rho     Density [kg/m^3]
%   h2.Cp      Specific heat [J/(kg.K)]
%   h2.k       Thermal conductivity [W/(m.K)]
%   h2.mu      Dynamic viscosity [Pa.s]
%   h2.nu      Kinematic viscosity [m^2/s]
%   h2.alpha   Thermal diffusivity [m^2/s]
%   h2.Pr      Prandtl number [-]
%   h2.beta    Thermal expansion coefficient [1/K]
%
% Valid for:
%   0 <= Tavg <= 200 °C
%
% Properties correspond to dry hydrogen at 1 atm.

%% Input checking

validateattributes(Tavg,{'numeric'},...
    {'real','finite','nonempty'},mfilename,'Tavg');

if any(Tavg(:)<0 | Tavg(:)>200)
    error('Tavg must be between 0 and 200 °C.');
end

%% Temperature

hydrogen.T_C = Tavg;
hydrogen.T_K = Tavg + 273.15;
hydrogen.p   = 101325;

%% Temperature table (°C)

Ttab = (0:10:200)';

%% Density [kg/m^3]

rho_tab = [...
0.0899;
0.0868;
0.0838;
0.0810;
0.0784;
0.0759;
0.0736;
0.0715;
0.0694;
0.0675;
0.0657;
0.0640;
0.0624;
0.0608;
0.0594;
0.0580;
0.0567;
0.0554;
0.0542;
0.0531;
0.0520];

%% Specific heat [J/(kg.K)]

Cp_tab = [...
14310;
14320;
14330;
14340;
14350;
14360;
14370;
14380;
14390;
14400;
14410;
14420;
14430;
14440;
14450;
14460;
14470;
14480;
14490;
14500;
14510];

%% Thermal conductivity [W/(m.K)]

k_tab = [...
0.168;
0.175;
0.182;
0.189;
0.196;
0.203;
0.210;
0.217;
0.224;
0.231;
0.238;
0.245;
0.252;
0.259;
0.266;
0.273;
0.280;
0.287;
0.294;
0.301;
0.308];

%% Dynamic viscosity [Pa.s]

mu_tab = [...
8.38e-6;
8.63e-6;
8.87e-6;
9.11e-6;
9.35e-6;
9.59e-6;
9.82e-6;
1.005e-5;
1.028e-5;
1.051e-5;
1.074e-5;
1.097e-5;
1.119e-5;
1.141e-5;
1.163e-5;
1.185e-5;
1.206e-5;
1.227e-5;
1.248e-5;
1.269e-5;
1.290e-5];

%% Prandtl number [-]

Pr_tab = [...
0.713;
0.711;
0.709;
0.707;
0.705;
0.703;
0.701;
0.699;
0.697;
0.695;
0.693;
0.691;
0.689;
0.687;
0.685;
0.683;
0.681;
0.679;
0.677;
0.675;
0.673];

%% Interpolation

hydrogen.rho = interp1(Ttab,rho_tab,Tavg,'pchip');
hydrogen.Cp  = interp1(Ttab,Cp_tab,Tavg,'pchip');
hydrogen.k   = interp1(Ttab,k_tab,Tavg,'pchip');
hydrogen.mu  = interp1(Ttab,mu_tab,Tavg,'pchip');
hydrogen.Pr  = interp1(Ttab,Pr_tab,Tavg,'pchip');

%% Derived properties

hydrogen.nu = hydrogen.mu ./ hydrogen.rho;

hydrogen.alpha = hydrogen.k ./ (hydrogen.rho .* hydrogen.Cp);

hydrogen.beta = 1 ./ hydrogen.T_K;

end
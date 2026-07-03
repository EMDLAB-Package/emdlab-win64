% EMDLAB: Electrical Machines Design Laboratory
% physical constants

classdef emdlab_phy_constants
    
    properties (Constant = true)
        
        % permeability of vaccume [H/m]
        mu0 = 4*pi*1e-7;  
        
        % reluctivity of vaccume [m/H]
        nu0 = 1/(4*pi*1e-7);  
        
        % permittivity of vaccume [H/m]
        e0 = 8.85e-12;  

        % Stefan-Boltzmann constant [W/(m^2 K^4)]
        sigmaSB = 5.670374419e-8;
        
    end
    
end
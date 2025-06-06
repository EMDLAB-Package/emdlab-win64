classdef emdlab_phy_material < handle
    
    properties
        
        % SI Unit
        %[W/(m.C)]
        ThermalConductivity (1,1) emdlab_phy_materialProperty;
        
        %[J/(Kg.C)]
        HeatCapacity (1,1) emdlab_phy_materialProperty;
        
        %[F/m]
        ElectricPermitivity (1,1) emdlab_phy_materialProperty;
        
        %[H/m]
        MagneticPermeability (1,1) emdlab_phy_materialProperty;
        
        ElectricConductivity (1,1) emdlab_phy_materialProperty;
        
        %[Kg/m^3]
        MassDensity (1,1) emdlab_phy_materialProperty;
        
        YoungModulus (1,1) emdlab_phy_materialProperty;
        
        PoissonRatio (1,1) emdlab_phy_materialProperty;       
        
    end
    
end

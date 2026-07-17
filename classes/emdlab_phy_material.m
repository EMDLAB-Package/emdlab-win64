% EMDLAB: Electrical Machines Design Laboratory
% all material properties must be in SI units

classdef emdlab_phy_material < handle

    properties

        % [W/(m.K)] or [W/(m.C)]
        ThermalConductivity (1,1) emdlab_phy_mp_ThermalConductivity;

        % [J/(Kg.C)]
        HeatCapacity (1,1) emdlab_phy_materialProperty;

        % [F/m]
        ElectricPermitivity (1,1) emdlab_phy_materialProperty;

        % [H/m]
        MagneticPermeability (1,1) emdlab_phy_materialProperty;

        % [S/m]
        ElectricConductivity (1,1) emdlab_phy_materialProperty;

        % [Kg/m^3]
        MassDensity (1,1) emdlab_phy_materialProperty;

        % [N/m^2] or [Pa]
        YoungModulus (1,1) emdlab_phy_materialProperty;

        % [-]
        PoissonRatio (1,1) emdlab_phy_materialProperty;

    end

    methods

        function setThermalConductivity(obj, value, unit)
            if nargin < 2, error('Not enough input arguments'); end
            if nargin < 3, unit = 'W/mK'; end
            if nargin > 3, error('Too many input arguments'); end
            obj.ThermalConductivity = obj.ThermalConductivity.setValue(value, unit);
        end

        function setMassDensity(obj, value)
            obj.MassDensity.value = value;
        end

        function setHeatCapacity(obj, value)
            obj.HeatCapacity.value = value;
        end

    end

end

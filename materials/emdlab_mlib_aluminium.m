classdef emdlab_mlib_aluminium < emdlab_phy_material

    methods

        function obj = emdlab_mlib_aluminium()

            % Aluminum material properties
            obj.ThermalConductivity.setValue(237,'W/mK');          % W/(m·K)
            obj.HeatCapacity.value = 900;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 3.5e7;       % S/m

            obj.MagneticPermeability.value = 1.0*4*pi*1e-7; % H/m (μ_r ≈ 1)
            obj.MagneticPermeability.isLinear = true;
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 2700;                 % kg/m³
            obj.YoungModulus.value = 7.0e10;              % Pa (≈ 70 GPa)
            obj.PoissonRatio.value = 0.33;                % dimensionless

        end

    end

end
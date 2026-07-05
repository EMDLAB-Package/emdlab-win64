classdef emdlab_mlib_copper < emdlab_phy_material

    methods

        function obj = emdlab_mlib_copper()

            % Copper material properties
            obj.ThermalConductivity.setValue(385,'W/mK');         % W/(m·K)
            obj.HeatCapacity.value = 385;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, copper ignored for quasi-static)
            obj.ElectricConductivity.value = 5.96e7;      % S/m

            obj.MagneticPermeability.value = 1.0*4*pi*1e-7;  % H/m (μ_r ≈ 1)
            obj.MagneticPermeability.isLinear = true;
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 8940;                 % kg/m³
            obj.YoungModulus.value = 1.1e11;              % Pa (≈ 110 GPa)
            obj.PoissonRatio.value = 0.34;                % dimensionless

        end

    end

end
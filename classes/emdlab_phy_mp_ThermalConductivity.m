% EMDLAB: Electrical Machines Design Laboratory
% emdlab -> physics -> material property -> thermal conductivity

classdef emdlab_phy_mp_ThermalConductivity

    properties(SetAccess = ?emdlab_phy_material)

        % value of material property -> thermal conductivity
        value (1,1) double;

        % this unit multiplier maps the value in original unit to SI unit
        unitUser (1,1) string; 
        unitSI (1,1) string = "W/(m.K)";
        unitMultiplier (1,1) double;

        % temperature dependent vs non-dependent
        % whe it is temperature dependent it should be defined as a function handle
        isTemperatureDependent (1,1) logical;

        % isotropic vs non-isotropic
        isIsotropic (1,1) logical;

        % homogeneous vs non-isotropic
        isHomogeneous (1,1) logical;

    end

    methods

        function obj = emdlab_phy_mp_ThermalConductivity()

            obj.value = 1;
            obj.isTemperatureDependent = false;
            obj.isIsotropic = true;
            obj.isHomogeneous = true;

        end

        function obj = setValue(obj, xValue, xUnit)

            % set default unit
            if nargin < 2, xUnit = 'W/(m.K)'; end

            obj.unitUser = erase(xUnit, ' ');
            xUnit = lower(obj.unitUser);
            % set unit multiplier
            switch xUnit

                case {'w/mk', 'w/(mk)', 'w/(m.k)', 'w/mc', 'w/(mc)', 'w/(m.c)'}
                    obj.unitMultiplier = 1;

                otherwise
                    error('Unsupported thermal conductivity unit: %s', obj.unitUser);

            end

            % set thermal conductivity
            if isnumeric(xValue) && isscalar(xValue)

                obj.value = xValue;
                obj.isTemperatureDependent = false;
                obj.isIsotropic = true;
                obj.isHomogeneous = true;

            elseif isnumeric(xValue) && isvector(xValue) && (length(xValue) == 3)

                obj.value = xValue(:); % xValue must be a [3x1] vector
                obj.isTemperatureDependent = false;
                obj.isIsotropic = false;
                obj.isHomogeneous = true;

            elseif isscalar(xValue) && isa(xValue,'function_handle')

                if nargin(xValue) == 1
                    obj.value = xValue;
                    obj.isTemperatureDependent = true;

                else
                end
                obj.isIsotropic = true;

            else
                error('Improper setting of thermal conductivity.');
            end

        end

        function y = getValue(obj, varargin)

            switch [obj.isTemperatureDependent, obj.isIsotropic, obj.isHomogeneous]

                case [false, true, true]
                    y = obj.value * obj.unitMultiplier;

                case [true, true, true]
                    if numel(varargin) == 1
                        if isscalar(varargin{1}) && isnumeric(varargin{1})
                            y = obj.value(varargin{1}) * obj.unitMultiplier;
                        else
                            error('Specified temperature must be a scalar numeric value.')
                        end
                    else
                        error(['When thermal conductivity is temperature dependent, ', ...
                            'to get thermal conductivity you must specify temperature.']);
                    end

                case [false, true, true]
                case [true, true, true]
                case [false, true, true]
                case [true, true, true]
                case [false, true, true]
                case [true, true, true]

            end

        end

    end
end

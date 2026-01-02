% EMDLAB: Electrical Machines Design Laboratory
% magnetic-static two-dimensional tl3
% triangular lagrangian elements common properties

classdef emdlab_solvers_ms2d_tlcp < handle

    properties (SetAccess = protected)

        % solver mesh
        m (1,1) emdlab_m2d_tmdb;

        % boundary conditions
        bcs (1,1) emdlab_bcs_scalarNodes;

        % elements data
        edata (1,1) struct;

        % results
        results (1,1) struct;

        % coils
        coils (1,1) struct;

        % coil arms
        coilArms (1,:) string;

    end

    properties (SetAccess = protected)

        % depth of problem
        depth (1,1) double {mustBePositive} = 1;

        % units
        units (1,1) emdlab_phy_units;

        % physical constants
        pcts (1,1) emdlab_phy_constants = emdlab_phy_constants;

        % solver Properties
        solverSettings (1,1) struct
        solverHistory (1,1) struct
        monitorResiduals (1,1) logical = false;
        solverIndex = 1;

        % states
        isBeEvaluated (1,1) logical = false;
        isBnEvaluated (1,1) logical = false;
        isElementDataAssigned (1,1) logical = false;
        isResultsValid (1,1) logical = false;

    end

    properties (Dependent = true)

        % number of coil arms
        NcoilArms (1,1) double;

        % number of coils
        Ncoils (1,1) double;

    end

    methods
        %% constructor and destructor
        function obj = emdlab_solvers_ms2d_tlcp()
            % default values
            obj.depth = 1;
            obj.bcs = emdlab_bcs_scalarNodes('TL3');
            obj.units = emdlab_phy_units;
        end

        function delete(obj)
            delete(obj.m);
            delete(obj.bcs);
            delete(obj.units);
        end

        function setLengthUnit(obj, unitValue)
            obj.units.setQuantityUnit('length', unitValue);
        end

        function setUnit(obj, varargin)
            obj.units.setQuantityUnit(varargin{:});
        end

        function setDepth(obj, value)
            obj.depth = value;
        end

        function setSolverIndex(obj, value)
            obj.solverIndex = value;
        end

        function y = get.NcoilArms(obj)
            y = numel(obj.coilArms);
        end

        function y = get.Ncoils(obj)
            y = numel(fieldnames(obj.coils));
        end

        %% solver properties for mesh zones
        % set default properties of a mesh zone
        function setdp(obj, mzName)

            obj.m.mzs.(mzName).props.isExcited = false;
            obj.m.mzs.(mzName).props.isCoilArm = false;
            obj.m.mzs.(mzName).props.isCoilMember = false;
            obj.m.mzs.(mzName).props.isMagnetized = false;
            obj.makeFalse_isElementDataAssigned;

        end

        function addmz(obj, mzName, mzValue)

            obj.m.addmz(mzName, mzValue);
            obj.setdp(mzName);
            obj.makeFalse_isElementDataAssigned;

        end

        function removemz(obj, mzName)

            obj.m.removemz(mzName);
            obj.makeFalse_isElementDataAssigned;

        end

        function rotateMeshZone(obj, mzName, varargin)

            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.m.rotateMeshZone(mzName, varargin{:});

            if obj.m.mzs.(mzName).props.isMagnetized
                obj.m.mzs.(mzName).props.magnetization.rotate(varargin{:});
            end

        end

        %% excitations definitions
        % returns coil name if coil exist otherwise it will raise an error
        function coilName = checkCoilExistence(obj, coilName)

            if ~isfield(obj.coils, coilName)
                throw(MException('', ['Coil with name [', coilName, '] does not exist.']));
            end

        end

        % returns coil name if coil does not exist otherwise it will raise an error
        function coilName = checkCoilNonExistence(obj, coilName)

            if isfield(obj.coils, coilName)
                throw(MException('', ['Another coil with name [', coilName, '] already exist.']));
            end

        end

        % define a new coil
        function defineCoil(obj, coilName)

            coilName = obj.checkCoilNonExistence(coilName);
            obj.coils.(coilName) = emdlab_solvers_ms2d_coil();
            obj.coils.(coilName).ci = obj.Ncoils;

        end

        % this function add a mesh zone to a coil and make it as a coil arm
        function addMeshZone2Coil(obj, coilName, mzName, turns, direction, kfill)

            % default arguments
            if nargin < 4, turns = 1; end
            if nargin < 5, direction = 1; end
            if nargin < 6, kfill = 1; end

            coilName = obj.checkCoilExistence(coilName);

            % loop over mesh zone names vector
            mzName = string(mzName);
            Nmzs = numel(mzName);

            if isscalar(turns)
                turns = turns*ones(1,Nmzs);
            else
                if length(turns) ~= Nmzs
                    error('Length of <turns> array must be the same as the number of mesh zones.');
                end
            end

            if isscalar(direction)
                direction = direction*ones(1,Nmzs);
            else
                if length(direction) ~= Nmzs
                    error('Length of <direction> array must be the same as the number of mesh zones.');
                end
            end

            if isscalar(kfill)
                kfill = kfill*ones(1,Nmzs);
            else
                if length(kfill) ~= Nmzs
                    error('Length of <kfill> array must be the same as the number of mesh zones.');
                end
            end

            for i = 1:Nmzs

                mzName(i) = obj.m.checkMeshZoneExistence(mzName(i));

                if obj.m.mzs.(mzName(i)).props.isCoilArm
                    error('Specified mesh zone is already defined as a coil arm.');
                end

                % get coil pointer
                cptr = obj.coils.(coilName);

                if ~ismember(direction(i), [-1,1])
                    error('The coil arm reference direction must be <1> or <-1>.');
                end

                if (kfill(i) < 0) || (kfill(i) > 1)
                    error('The coil fill factor must between zero to one.');
                end

                cptr.addCoilArm(mzName(i), direction(i));
                obj.m.mzs.(mzName(i)).props.turns = turns(i);
                obj.m.mzs.(mzName(i)).props.direction = direction(i);
                obj.m.mzs.(mzName(i)).props.kfill = kfill(i);
                obj.m.mzs.(mzName(i)).props.isCoilArm = true;
                obj.coilArms(end+1) = mzName(i);
                obj.m.mzs.(mzName(i)).props.cai = obj.NcoilArms;

            end

        end

        % this function add a number of mesh zones to a coil and make it as a coil arm
        function addMeshZones2Coil(obj, coilName, mzName, turns, kfill)

            % default arguments
            if nargin < 4, turns = 1; end
            if nargin < 5, kfill = 1; end

            coilName = obj.checkCoilExistence(coilName);

            % loop over mesh zone names vector
            mzName = string(mzName);
            Nmzs = numel(mzName);

            if isscalar(turns)
                turns = turns*ones(1,Nmzs);
            else
                if length(turns) ~= Nmzs
                    error('Length of <turns> array must be the same as the number of mesh zones.');
                end
            end

            if isscalar(kfill)
                kfill = kfill*ones(1,Nmzs);
            else
                if length(kfill) ~= Nmzs
                    error('Length of <kfill> array must be the same as the number of mesh zones.');
                end
            end

            for i = 1:Nmzs

                mzName(i) = obj.m.checkMeshZoneExistence(mzName(i));

                if obj.m.mzs.(mzName(i)).props.isCoilArm
                    error('Specified mesh zone is already defined as a coil arm.');
                end

                % get coil pointer
                cptr = obj.coils.(coilName);

                direction = 1;
                if turns(i)<0, direction = -1; end

                if (kfill(i) < 0) || (kfill(i) > 1)
                    error('The coil fill factor must between zero to one.');
                end

                cptr.addCoilArm(mzName(i), direction);
                obj.m.mzs.(mzName(i)).props.turns = abs(turns(i));
                obj.m.mzs.(mzName(i)).props.direction = direction;
                obj.m.mzs.(mzName(i)).props.kfill = kfill(i);
                obj.m.mzs.(mzName(i)).props.isCoilArm = true;
                obj.coilArms(end+1) = mzName(i);
                obj.m.mzs.(mzName(i)).props.cai = obj.NcoilArms;

            end

        end

        % check if all defined coils have coil arms, and also calculates Rdc
        function checkCoils(obj)

            % coils with zeros number of coil arms
            coilNames = fieldnames(obj.coils);
            for i = 1:obj.Ncoils

                if obj.coils.(coilNames{i}).NcoilArms == 0
                    error(['Coil <', coilNames{i}, '> does not have any coil arm.']);
                end

            end

        end

        % setting coil current
        function setCoilCurrent(obj, coilName, value)
            coilName = obj.checkCoilExistence(coilName);
            obj.coils.(coilName).setCurrent(value);
            if obj.isElementDataAssigned
                % applying current of excitation matrices
                coilNames = fieldnames(obj.coils);

                for i = 1:numel(coilNames)

                    % get coil pointer
                    cptr = obj.coils.(coilNames{i});

                    for j = 1:cptr.NcoilArms

                        % get coil arm pointer
                        mzptr = obj.m.mzs.(cptr.coilArms(j));

                        % set coil arm current density
                        obj.edata.InternalCurrentDensity(:,obj.m.ezi(:, mzptr.zi)) = ...
                            mzptr.props.direction * mzptr.props.turns * cptr.current * obj.units.k_current ...
                            / (mzptr.getArea * obj.units.k_length^2);

                    end

                end
            end

        end

        % definition of excitation that can be current or current density
        function setExcitation(obj, mzName, varargin)

            mzName = obj.m.checkMeshZoneExistence(mzName);

            if obj.m.mzs.(mzName).props.isWindingMember
                throw(MException('', ['Mesh zone [', mzName, '] already is assinged to a winding.']));
            end

            obj.m.mzs.(mzName).props.excitation = emdlab_solvers_ms2d_excitation(varargin{:});
            obj.m.mzs.(mzName).props.isExcited = true;
            % change states
            obj.makeFalse_isElementDataAssigned;

        end

        % this function set the magnetization property of a mesh zone
        function setMagnetization(obj, mzName, varargin)
            mzName = obj.m.checkMeshZoneExistence(mzName);

            if nargin == 3 && isa(varargin{1}, 'emdlab_solvers_ms2d_magnetization')
                obj.m.mzs.(mzName).props.magnetization = varargin{1};
            else
                obj.m.mzs.(mzName).props.magnetization = emdlab_solvers_ms2d_magnetization(varargin{:});
            end

            obj.m.mzs.(mzName).props.isMagnetized = true;
            % change states
            obj.makeFalse_isElementDataAssigned;
        end

        %% visualization functions
        % show a single coils
        function varargout = showCoil(obj, coilName, varargin)

            switch numel(varargin)
                case 0
                    f = figure;
                    ax = axes(f);
                    f.Name = 'Coil arm';

                case 1

                    if isa(varargin{1}, 'matlab.graphics.axis.Axes')
                        ax = varargin{1};
                        f = ax.Parent;
                    end

                otherwise
                    error('Too many input arguments');
            end

            coilName = obj.checkCoilExistence(coilName);
            cptr = obj.coils.(coilName);

            for i = 1:cptr.NcoilArms
                mzptr = obj.m.mzs.(cptr.coilArms(i));

                if cptr.directions(i) == 1
                    mzColor = 'b';
                else
                    mzColor = 'r';
                end

                if mzptr.props.turns == 0
                    mzColor = [0.8,0.8,0.8];
                end

                patch('Faces', mzptr.cl, 'Vertices', mzptr.nodes, ...
                    'FaceColor', mzColor, 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'parent', ax);
            end

            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        % show all coil arms
        function varargout = showAllCoilArms(obj, varargin)

            switch numel(varargin)
                case 0
                    f = figure;
                    ax = axes(f);
                    f.Name = 'All defined coil arms';

                case 1

                    if isa(varargin{1}, 'matlab.graphics.axis.Axes')
                        ax = varargin{1};
                        f = ax.Parent;
                    end

                otherwise
                    error('Too many input arguments');
            end

            mzNames = fieldnames(obj.m.mzs);

            for i = 1:numel(mzNames)
                mzptr = obj.m.mzs.(mzNames{i});

                if mzptr.props.isCoilArm
                    patch('Faces', mzptr.cl, 'Vertices', mzptr.nodes, ...
                        'FaceColor', 'g', 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'parent', ax);
                end

            end

            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        %% boundary conditions
        function clearAllBCs(obj)
            obj.bcs.clearAllBCs;
        end

        function setAzBC(obj, index, value, varargin)
            obj.bcs.setDirichlet(index, value, varargin{:});
        end

        function setOddPeriodicBC(obj, varargin)
            obj.bcs.setOddPeriodic(varargin{:});
        end

        function setEvenPeriodicBC(obj, varargin)
            obj.bcs.setEvenPeriodic(varargin{:});
        end

        %% post-proccessing: field quantities calculations
        function clearAllResults(obj)

            % result names
            rNames = fieldnames(obj.results);

            for i = 1:numel(rNames)
                obj.results.(rNames{i}) = [];
            end

            obj.isResultsValid = false;

        end

        % Evaluation of B [tesla] on gaussian and mesh points of each element
        function evalBe(obj)

            if obj.m.isTL3

                % 1 gaussian point -> (1/3,1/3)
                % 3 mesh points at corners
                [obj.results.Bxg, obj.results.Byg, obj.results.Bxn, obj.results.Byn] = ...
                    emdlab_m2d_tl3_evalB(obj.m.cl, obj.results.A, obj.m.JIT);
                obj.results.Bxg = obj.results.Bxg * (obj.units.k_magneticVectorPotential / obj.units.k_length);
                obj.results.Byg = obj.results.Byg * (obj.units.k_magneticVectorPotential / obj.units.k_length);
                obj.results.Bxn = obj.results.Bxn * (obj.units.k_magneticVectorPotential / obj.units.k_length);
                obj.results.Byn = obj.results.Byn * (obj.units.k_magneticVectorPotential / obj.units.k_length);

            elseif obj.m.isTL6

                % Evaluation of B on gaussian points
                [obj.results.Bxg, obj.results.Byg, obj.results.Bxn, obj.results.Byn] = ...
                    emdlab_m2d_tl6_evalB(obj.m.cl, obj.results.A, obj.m.JIT);
                obj.results.Bxg = obj.results.Bxg * (obj.units.k_magneticVectorPotential / obj.units.k_length);
                obj.results.Byg = obj.results.Byg * (obj.units.k_magneticVectorPotential / obj.units.k_length);
                obj.results.Bxn = obj.results.Bxn * (obj.units.k_magneticVectorPotential / obj.units.k_length);
                obj.results.Byn = obj.results.Byn * (obj.units.k_magneticVectorPotential / obj.units.k_length);

            else
                error('Wrong element type');
            end

            obj.isBeEvaluated = true;

        end

        % Evaluation of H [A/m] on gaussian and mesh points of each element
        function evalHe(obj)

            obj.results.Hxg = obj.edata.MagneticReluctivity .* obj.results.Bxg - obj.edata.MagnetizationX;
            obj.results.Hyg = obj.edata.MagneticReluctivity .* obj.results.Byg - obj.edata.MagnetizationY;

            obj.results.Hxn = obj.edata.MagneticReluctivity .* obj.results.Bxn - obj.edata.MagnetizationX;
            obj.results.Hyn = obj.edata.MagneticReluctivity .* obj.results.Byn - obj.edata.MagnetizationY;

        end

        % evaluate smoothed B on mesh points
        function evalBn(obj)

            obj.evalBe;

            mzsNames = fieldnames(obj.m.mzs);
            obj.results.BxnSmooth = zeros(3,obj.m.Ne);
            obj.results.BynSmooth = zeros(3,obj.m.Ne);

            for i = 1:numel(mzsNames)

                mzptr = obj.m.mzs.(mzsNames{i});
                eziptr = obj.m.ezi(:,mzptr.zi);

                if obj.m.isTL3

                    [obj.results.BxnSmooth(:,eziptr), obj.results.BynSmooth(:,eziptr)] = ...
                        emdlab_m2d_tl3_evalBnSmooth(obj.m.cl(eziptr,:), obj.results.Bxn(:,eziptr), obj.results.Byn(:,eziptr), obj.m.gea(eziptr), obj.m.Nn);

                elseif obj.m.isTL6

                    [obj.results.BxnSmooth(:,eziptr), obj.results.BynSmooth(:,eziptr)] = ...
                        emdlab_m2d_tl6_evalBnSmooth(obj.m.cl(eziptr,:), obj.results.Bxn(:,eziptr), obj.results.Byn(:,eziptr), obj.m.gea(eziptr), obj.m.Nn);

                else
                    error('Wrong element type');
                end

            end

            obj.isBnEvaluated = true;

        end

        % evaluate smoothed H on mesh points
        function evalHn(obj)

            obj.evalHe;

            mzsNames = fieldnames(obj.m.mzs);
            obj.results.HxnSmooth = zeros(3,obj.m.Ne);
            obj.results.HynSmooth = zeros(3,obj.m.Ne);

            for i = 1:numel(mzsNames)

                mzptr = obj.m.mzs.(mzsNames{i});
                eziptr = obj.m.ezi(:,mzptr.zi);

                if obj.m.isTL3

                    [obj.results.HxnSmooth(:,eziptr), obj.results.HynSmooth(:,eziptr)] = ...
                        emdlab_m2d_tl3_evalBnSmooth(obj.m.cl(eziptr,:), obj.results.Hxn(:,eziptr), obj.results.Hyn(:,eziptr), obj.m.gea(eziptr), obj.m.Nn);

                elseif obj.m.isTL6

                    [obj.results.HxnSmooth(:,eziptr), obj.results.HynSmooth(:,eziptr)] = ...
                        emdlab_m2d_tl6_evalBnSmooth(obj.m.cl(eziptr,:), obj.results.Hxn(:,eziptr), obj.results.Hyn(:,eziptr), obj.m.gea(eziptr), obj.m.Nn);

                else
                    error('Wrong element type');
                end


            end

        end

        %% post-proccessing: energy calculations
        % calculate total stored energy and co-energy
        function [ye, yc] = evalTotalEnergyCoenergyTL3(obj)

            mzsNames = fieldnames(obj.m.mzs);
            Bk2 = obj.results.Bxg.^2 + obj.results.Byg.^2;
            ye = 0;
            yc = 0;

            for i = 1:numel(mzsNames)

                mzptr = obj.m.mzs.(mzsNames{i});
                mptr = obj.m.mts.(mzptr.material);
                eziptr = obj.m.ezi(:,mzptr.zi);

                if mptr.MagneticPermeability.isLinear

                    % mesh zone energy
                    mze = 0.5*obj.edata.MagneticReluctivity(eziptr) * (mzptr.getAreaOfElements .* Bk2(eziptr)');
                    ye = ye + mze;
                    yc = yc + mze;

                else

                    % mesh zone energy
                    Bk = sqrt(Bk2(eziptr)');
                    ye = ye + mzptr.getAreaOfElements' * ppval(mptr.weB, Bk);
                    yc = yc + mzptr.getAreaOfElements' * ppval(mptr.wcH, ppval(mptr.HB, Bk));

                end

            end

            ye = ye*obj.getDepth*obj.units.k_length^2;
            yc = yc*obj.getDepth*obj.units.k_length^2;

        end

        function [ye, yc] = evalTotalEnergyCoenergyTL6(obj)

            mzsNames = fieldnames(obj.m.mzs);
            Bk2 = obj.results.Bxg.^2 + obj.results.Byg.^2;
            ye = 0;
            yc = 0;

            for i = 1:numel(mzsNames)

                mzptr = obj.m.mzs.(mzsNames{i});
                mptr = obj.m.mts.(mzptr.material);
                eziptr = obj.m.ezi(:,mzptr.zi);

                if mptr.MagneticPermeability.isLinear

                    % mesh zone energy
                    mze = mzptr.getAreaOfElements' .* (obj.edata.MagneticReluctivity(:,eziptr) .* Bk2(:,eziptr));
                    mze = sum(sum(mze))/6;
                    ye = ye + mze;
                    yc = yc + mze;

                else

                    % mesh zone energy
                    Bk = sqrt(Bk2(:,eziptr));
                    mze = mzptr.getAreaOfElements' .* ppval(mptr.weB, Bk);
                    mze = sum(sum(mze))/3;
                    ye = ye + mze;
                    mzc = mzptr.getAreaOfElements' .* ppval(mptr.wcH, ppval(mptr.HB, Bk));
                    mzc= sum(sum(mzc))/3;
                    yc = yc + mzc;

                end

            end

            ye = ye*obj.getDepth*obj.units.k_length^2;
            yc = yc*obj.getDepth*obj.units.k_length^2;

        end

        function [ye, yc] = evalTotalEnergyCoenergy(obj)

            if obj.m.isTL3
                [ye, yc] = obj.evalTotalEnergyCoenergyTL3;
            elseif obj.m.isTL6
                [ye, yc] = obj.evalTotalEnergyCoenergyTL6;
            else
                error('Wrong element type');
            end

        end

        %% post-proccessing: flux linkage and inductance matrix calculations
        % calculate mesh zone flux linkage
        function y = evalMeshZoneFluxLinkageTL3(obj, mzName)

            mzName = obj.m.checkMeshZoneExistence(mzName);
            y = obj.m.mzs.(mzName).getQ * obj.results.A(obj.m.mzs.(mzName).l2g) * obj.getDepth;

        end

        function y = evalMeshZoneFluxLinkageTL6(obj, mzName)

            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            eziptr = obj.m.ezi(:,mzptr.zi);

            y = obj.m.gea(eziptr) * sum(obj.results.A(obj.m.cl(eziptr,4:6)),2) * (obj.getDepth/3) / mzptr.getArea;

        end

        function y = evalMeshZoneFluxLinkage(obj, mzName)

            if obj.m.isTL3
                y = obj.evalMeshZoneFluxLinkageTL3(mzName);
            elseif obj.m.isTL6
                y = obj.evalMeshZoneFluxLinkageTL6(mzName);
            else
                error('Wrong element type');
            end

        end

        % calculate coil flux linkage
        function y = evalCoilFluxLinkage(obj, coilName)

            coilName = obj.checkCoilExistence(coilName);
            cptr = obj.coils.(coilName);
            y = 0;

            for i = 1:cptr.NcoilArms
                % get pointer to coil arm
                mzptr = obj.m.mzs.(cptr.coilArms(i));

                % eval and add mesh zone flux linkage
                y = y + mzptr.props.direction * mzptr.props.turns * obj.evalMeshZoneFluxLinkage(cptr.coilArms(i));
            end

            cptr.fluxLinkage = y;

        end

        %% post-proccessing: torque and force calculations
        % torque evaluation using maxwell stress tensor
        function y = evalTorqueByMST(obj, xc, yc, r, N)

            if nargin < 5, N = 10001; end
            [br, bt, t] = obj.getBrBtOnCircle(xc, yc, r, N);
            br(isnan(br)) = 0;
            bt(isnan(bt)) = 0;
            y = trapz(t, br.*bt) * r^2 * obj.units.k_length^2 * obj.getDepth / (4 * pi * 1e-7);

        end

        function [Ft,Fn] = evalForceByMSTOnLine(obj, x1, y1, x2, y2, N)

            if nargin < 6, N = 10001; end
            x = linspace(x1,x2,N);
            y = linspace(y1,y2,N);
            [Bx, By] = obj.getBxByOnPoints(x,y);
            Bx(isnan(Bx)) = 0;
            By(isnan(By)) = 0;
            t = [x2-x1, y2-y1];
            t = t/norm(t);
            n = [-t(2),t(1)];
            Bt = Bx*t(1) + By*t(2);
            Bn = Bx*n(1) + By*n(2);
            l = cumsum([0,sqrt(diff(x).^2+diff(y).^2)]);
            Ft = trapz(l, Bt.*Bn) * obj.units.k_length^2 * obj.getDepth / (4 * pi * 1e-7);
            Fn = trapz(l, 0.5*(Bn.^2 - Bt.^2)) * obj.units.k_length^2 * obj.getDepth / (4 * pi * 1e-7);

        end

        % torque evaluation using maxwell stress tensor on 3 layers
        function y = evalTorqueByMST3(obj, xc, yc, r, gap, N)

            if nargin < 6, N = 10001; end
            y = (obj.evalTorqueByMST(xc,yc,r-gap/4,N) + obj.evalTorqueByMST(xc,yc,r,N) + obj.evalTorqueByMST(xc,yc,r+gap/4,N))/3;

        end

        % torque evaluation using surface MST on a number of selected objects
        function torque = evalTorqueBySurfaceMSTTL3(obj, xc, yc, varargin)

            % get names string list
            mzNames = emdlab_flib_varargin2StringList(varargin{:});

            % index of selected zones
            zi = zeros(1,numel(mzNames));

            % check mesh zone existance
            for i = 1:numel(mzNames)
                mzNames(i) = obj.m.checkMeshZoneExistence(mzNames(i));
                zi(i) = obj.m.mzs.(mzNames(i)).zi;
            end

            % find boundary edges of mesh zone exposing to air
            eIndices = ismember(obj.m.edges(:, 3),zi) & (~ismember(obj.m.edges(:, 4),zi));
            eIndices = eIndices | (ismember(obj.m.edges(:, 4),zi) & (~ismember(obj.m.edges(:, 3),zi)));
            eIndices = eIndices & (~ obj.m.bedges);
            eIndices = find(eIndices);

            torque = 0;
            for eIndex = eIndices'

                p1Index = obj.m.edges(eIndex,1);
                p2Index = obj.m.edges(eIndex,2);

                r1 = obj.m.nodes(p1Index,:) - xc;
                r2 = obj.m.nodes(p2Index,:) - yc;

                % finding normal vector
                n = r2 - r1;
                el = norm(n); % edge length
                n = n/el;

                if ismember(obj.m.edges(eIndex, 3),zi)
                    n = ext_protate2(n, -pi/2);
                    elIndex = obj.m.edges(eIndex, 7);
                elseif ismember(obj.m.edges(eIndex, 4),zi)
                    n = ext_protate2(n, pi/2);
                    elIndex = obj.m.edges(eIndex, 5);
                else
                    error('Internal error.');
                end

                index1 = find(p1Index == obj.m.cl(elIndex,:));
                index2 = find(p2Index == obj.m.cl(elIndex,:));

                mu0 = 4*pi*1e-7;

                Bx1 = obj.results.BxnSmooth(index1, elIndex);
                By1 = obj.results.BynSmooth(index1, elIndex);

                Bx2 = obj.results.BxnSmooth(index2, elIndex);
                By2 = obj.results.BynSmooth(index2, elIndex);

                Txx = (0.5/mu0) * (Bx1^2 - By1^2);
                Txy_yx = (1/mu0) * (Bx1 * By1);
                Tyy = (0.5/mu0) * (By1^2 - Bx1^2);

                Fx1 = (Txx * n(1) + Txy_yx * n(2)) * el;
                Fy1 = (Txy_yx * n(1) + Tyy * n(2)) * el;

                Txx = (0.5/mu0) * (Bx2^2 - By2^2);
                Txy_yx = (1/mu0) * (Bx2 * By2);
                Tyy = (0.5/mu0) * (By2^2 - Bx2^2);

                Fx2 = (Txx * n(1) + Txy_yx * n(2)) * el;
                Fy2 = (Txy_yx * n(1) + Tyy * n(2)) * el;

                r21 = r2 - r1;
                F21 = [Fx2, Fy2] - [Fx1, Fy1];

                a = r1(1) * Fy1 - r1(2) * Fx1;
                b = r1(1) * F21(2) - r1(2) * F21(1) + r21(1) * Fy1 - r21(2) * Fx1;
                c = r21(1) * F21(2) - r21(2) * F21(1);

                torque = torque + a + b/2 + c/3;

                %                 Bx = obj.results.Bxg(elIndex);
                %                 By = obj.results.Byg(elIndex);
                %
                %                 mu0 = 4*pi*1e-7;
                %                 Txx = (0.5/mu0) * (Bx^2 - By^2);
                %                 Txy_yx = (1/mu0) * (Bx * By);
                %                 Tyy = (0.5/mu0) * (By^2 - Bx^2);
                %
                %                 Fx = (Txx * n(1) + Txy_yx * n(2)) * el;
                %                 Fy = (Txy_yx * n(1) + Tyy * n(2)) * el;
                %
                %                 a = r1(1) * Fy - r1(2) * Fx;
                %                 r = r2 - r1;
                %                 b = r(1) * Fy - r(2) * Fx;
                %
                %                 torque = torque + a + b/2;

            end

            torque = torque * obj.getDepth * obj.units.k_length^2;

        end

        function torque = evalTorqueBySurfaceMSTTL6(obj, xc, yc, varargin)

            % get names string list
            mzNames = emdlab_flib_varargin2StringList(varargin{:});

            % index of selected zones
            zi = zeros(1,numel(mzNames));

            % check mesh zone existance
            for i = 1:numel(mzNames)
                mzNames(i) = obj.m.checkMeshZoneExistence(mzNames(i));
                zi(i) = obj.m.mzs.(mzNames(i)).zi;
            end

            % find boundary edges of mesh zone exposing to air
            eIndices = ismember(obj.m.edges(:, 3),zi) & (~ismember(obj.m.edges(:, 4),zi));
            eIndices = eIndices | (ismember(obj.m.edges(:, 4),zi) & (~ismember(obj.m.edges(:, 3),zi)));
            eIndices = eIndices & (~ obj.m.bedges);
            eIndices = find(eIndices);

            torque = 0;
            for eIndex = eIndices'

                p1Index = obj.m.edges(eIndex,1);
                p2Index = obj.m.edges(eIndex,2);

                r1 = obj.m.nodes(p1Index,:) - xc;
                r2 = obj.m.nodes(p2Index,:) - yc;

                % finding normal vector
                n = r2 - r1;
                el = norm(n); % edge length
                n = n/el;

                if ismember(obj.m.edges(eIndex, 3),zi)
                    n = ext_protate2(n, -pi/2);
                    elIndex = obj.m.edges(eIndex, 7);
                elseif ismember(obj.m.edges(eIndex, 4),zi)
                    n = ext_protate2(n, pi/2);
                    elIndex = obj.m.edges(eIndex, 5);
                else
                    error('Internal error.');
                end

                index1 = find(p1Index == obj.m.cl(elIndex,:));
                index2 = find(p2Index == obj.m.cl(elIndex,:));

                mu0 = 4*pi*1e-7;

                Bx1 = obj.results.BxnSmooth(index1, elIndex);
                By1 = obj.results.BynSmooth(index1, elIndex);

                Bx2 = obj.results.BxnSmooth(index2, elIndex);
                By2 = obj.results.BynSmooth(index2, elIndex);

                Txx = (0.5/mu0) * (Bx1^2 - By1^2);
                Txy_yx = (1/mu0) * (Bx1 * By1);
                Tyy = (0.5/mu0) * (By1^2 - Bx1^2);

                Fx1 = (Txx * n(1) + Txy_yx * n(2)) * el;
                Fy1 = (Txy_yx * n(1) + Tyy * n(2)) * el;

                Txx = (0.5/mu0) * (Bx2^2 - By2^2);
                Txy_yx = (1/mu0) * (Bx2 * By2);
                Tyy = (0.5/mu0) * (By2^2 - Bx2^2);

                Fx2 = (Txx * n(1) + Txy_yx * n(2)) * el;
                Fy2 = (Txy_yx * n(1) + Tyy * n(2)) * el;

                r21 = r2 - r1;
                F21 = [Fx2, Fy2] - [Fx1, Fy1];

                a = r1(1) * Fy1 - r1(2) * Fx1;
                b = r1(1) * F21(2) - r1(2) * F21(1) + r21(1) * Fy1 - r21(2) * Fx1;
                c = r21(1) * F21(2) - r21(2) * F21(1);

                torque = torque + a + b/2 + c/3;

            end

            torque = torque * obj.getDepth * obj.units.k_length^2;

        end

        function torque = evalTorqueBySurfaceMST(obj, xc, yc, varargin)

            if obj.m.isTL3
                torque = obj.evalTorqueBySurfaceMSTTL3(xc, yc, varargin{:});
            elseif obj.m.isTL6
                torque = obj.evalTorqueBySurfaceMSTTL6(xc, yc, varargin{:});
            else
                error('Wrong element type');
            end

        end

        function [x,y,Fx,Fy] = evalSurfaceForceDensityMST(obj, varargin)

            % get names string list
            mzNames = emdlab_flib_varargin2StringList(varargin{:});

            % index of selected zones
            zi = zeros(1,numel(mzNames));

            % check mesh zone existance
            for i = 1:numel(mzNames)
                mzNames(i) = obj.m.checkMeshZoneExistence(mzNames(i));
                zi(i) = obj.m.mzs.(mzNames(i)).zi;
            end

            % find boundary edges of mesh zone exposing to air
            eIndices = ismember(obj.m.edges(:, 3),zi) & (~ismember(obj.m.edges(:, 4),zi));
            eIndices = eIndices | (ismember(obj.m.edges(:, 4),zi) & (~ismember(obj.m.edges(:, 3),zi)));
            eIndices = eIndices & (~ obj.m.bedges);
            eIndices = find(eIndices);

            Ne = length(eIndices);
            x = zeros(1,Ne);
            y = zeros(1,Ne);
            Fx = zeros(1,Ne);
            Fy = zeros(1,Ne);
            cnt = 0;

            torque = 0;
            for eIndex = eIndices'

                p1Index = obj.m.edges(eIndex,1);
                p2Index = obj.m.edges(eIndex,2);

                r1 = obj.m.nodes(p1Index,:);
                r2 = obj.m.nodes(p2Index,:);

                % finding normal vector
                n = r2 - r1;
                el = norm(n); % edge length
                n = n/el;

                if ismember(obj.m.edges(eIndex, 3),zi)
                    n = ext_protate2(n, -pi/2);
                    elIndex = obj.m.edges(eIndex, 7);
                elseif ismember(obj.m.edges(eIndex, 4),zi)
                    n = ext_protate2(n, pi/2);
                    elIndex = obj.m.edges(eIndex, 5);
                else
                    error('Internal error.');
                end

                index1 = find(p1Index == obj.m.cl(elIndex,:));
                index2 = find(p2Index == obj.m.cl(elIndex,:));

                mu0 = 4*pi*1e-7;

                Bx1 = obj.results.BxnSmooth(index1, elIndex);
                By1 = obj.results.BynSmooth(index1, elIndex);

                Bx2 = obj.results.BxnSmooth(index2, elIndex);
                By2 = obj.results.BynSmooth(index2, elIndex);

                Txx = (0.5/mu0) * (Bx1^2 - By1^2);
                Txy_yx = (1/mu0) * (Bx1 * By1);
                Tyy = (0.5/mu0) * (By1^2 - Bx1^2);

                Fx1 = (Txx * n(1) + Txy_yx * n(2)) * el;
                Fy1 = (Txy_yx * n(1) + Tyy * n(2)) * el;

                Txx = (0.5/mu0) * (Bx2^2 - By2^2);
                Txy_yx = (1/mu0) * (Bx2 * By2);
                Tyy = (0.5/mu0) * (By2^2 - Bx2^2);

                Fx2 = (Txx * n(1) + Txy_yx * n(2)) * el;
                Fy2 = (Txy_yx * n(1) + Tyy * n(2)) * el;

                r21 = r2 - r1;
                F21 = [Fx2, Fy2] - [Fx1, Fy1];

                a = r1(1) * Fy1 - r1(2) * Fx1;
                b = r1(1) * F21(2) - r1(2) * F21(1) + r21(1) * Fy1 - r21(2) * Fx1;
                c = r21(1) * F21(2) - r21(2) * F21(1);

                torque = torque + a + b/2 + c/3;

                rmid = (r1+r2)/2;

                cnt = cnt + 1;
                x(cnt) = rmid(1);
                y(cnt) = rmid(2);
                Fx(cnt) = (Fx1+Fx2)/2;
                Fy(cnt) = (Fy1+Fy2)/2;


            end

        end

        % torque evaluation by Arrkio method
        function torque = evalTorqueByArkkioTL3(obj, mzName, GapLength)

            obj.evalBe;
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            b = [obj.results.Bxg(obj.m.ezi(:, mzptr.zi)); obj.results.Byg(obj.m.ezi(:, mzptr.zi))]';
            p = mzptr.getCenterOfElements;
            b = ext_xy2rt(p, b);
            r = sqrt(sum(p.^2, 2));
            torque = mzptr.getAreaOfElements' * (b(:, 1) .* b(:, 2) .* r);
            torque = torque * obj.getDepth * obj.units.k_length^2 / GapLength / (4 * pi * 1e-7);

        end

        function torque = evalTorqueByArkkioTL6(obj, mzName, GapLength)

            obj.evalBe;
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);

            % barycentric coordinates
            k1 = [0.5,0.5,0];
            k2 = [0.5,0,0.5];
            k3 = [0,0.5,0.5];

            % weights
            w = [1/3,1/3,1/3];

            % x and y coordinates of mesh zone node elements
            cl = obj.m.cl(obj.m.ezi(:, mzptr.zi),[1,2,3]);
            x = obj.m.nodes(:,1); x = x(cl);
            y = obj.m.nodes(:,2); y = y(cl);

            torque = 0;
            for i = 1:3

                b = [obj.results.Bxg(i,obj.m.ezi(:, mzptr.zi)); obj.results.Byg(i,obj.m.ezi(:, mzptr.zi))]';

                p = [k1(i)*x(:,1) + k2(i)*x(:,2) + k3(i)*x(:,3), ...
                    k1(i)*y(:,1) + k2(i)*y(:,2) + k3(i)*y(:,3)];

                b = ext_xy2rt(p, b);
                r = sqrt(sum(p.^2, 2));
                torque_i = mzptr.getAreaOfElements' * (b(:, 1) .* b(:, 2) .* r)*w(i);
                torque = torque + torque_i * obj.getDepth * obj.units.k_length^2 / GapLength / (4 * pi * 1e-7);

            end

        end

        function torque = evalTorqueByArkkio(obj, mzName, GapLength)

            if obj.m.isTL3
                torque = obj.evalTorqueByArkkioTL3(mzName, GapLength);
            elseif obj.m.isTL6
                torque = obj.evalTorqueByArkkioTL6(mzName, GapLength);
            else
                error('Wrong element type');
            end

        end

        % torque evaluation by Arrkio method using smoothed field data
        function torque = evalTorqueByArkkioSmoothTL3(obj, mzName, GapLength)

            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            bx = sum(obj.results.BxnSmooth(:,obj.m.ezi(:, mzptr.zi)),1)'/3;
            by = sum(obj.results.BynSmooth(:,obj.m.ezi(:, mzptr.zi)),1)'/3;
            b = [bx,by];
            p = mzptr.getCenterOfElements;
            b = ext_xy2rt(p, b);
            r = sqrt(sum(p.^2, 2));
            torque = mzptr.getAreaOfElements' * (b(:, 1) .* b(:, 2) .* r);
            torque = torque * obj.getDepth * obj.units.k_length^2 / GapLength / (4 * pi * 1e-7);

        end

        function torque = evalTorqueByArkkioSmooth(obj, mzName, GapLength)

            if obj.m.isTL3
                torque = obj.evalTorqueByArkkioSmoothTL3(mzName, GapLength);
            elseif obj.m.isTL6
                torque = obj.evalTorqueByArkkioTL6(mzName, GapLength);
            else
                error('Wrong element type');
            end

        end

        % force evaluation by Arrkio method
        function torque = evalForceByArkkioTL3(obj, mzName, GapLength)

            obj.evalBe;
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            b = [obj.results.Bxg(obj.m.ezi(:, mzptr.zi)); obj.results.Byg(obj.m.ezi(:, mzptr.zi))]';
            p = mzptr.getCenterOfElements;
            b = ext_xy2rt(p, b);
            r = sqrt(sum(p.^2, 2));
            torque = mzptr.getAreaOfElements' * (b(:, 1) .* b(:, 2) .* r);
            torque = torque * obj.getDepth * obj.units.k_length^2 / GapLength / (4 * pi * 1e-7);

        end

        function torque = evalForceByArkkioTL6(obj, mzName, GapLength)

            obj.evalBe;
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);

            % barycentric coordinates
            k1 = [0.5,0.5,0];
            k2 = [0.5,0,0.5];
            k3 = [0,0.5,0.5];

            % weights
            w = [1/3,1/3,1/3];

            % x and y coordinates of mesh zone node elements
            cl = obj.m.cl(obj.m.ezi(:, mzptr.zi),[1,2,3]);
            x = obj.m.nodes(:,1); x = x(cl);
            y = obj.m.nodes(:,2); y = y(cl);

            torque = 0;
            for i = 1:3

                b = [obj.results.Bxg(i,obj.m.ezi(:, mzptr.zi)); obj.results.Byg(i,obj.m.ezi(:, mzptr.zi))]';

                p = [k1(i)*x(:,1) + k2(i)*x(:,2) + k3(i)*x(:,3), ...
                    k1(i)*y(:,1) + k2(i)*y(:,2) + k3(i)*y(:,3)];

                b = ext_xy2rt(p, b);
                r = sqrt(sum(p.^2, 2));
                torque_i = mzptr.getAreaOfElements' * (b(:, 1) .* b(:, 2) .* r)*w(i);
                torque = torque + torque_i * obj.getDepth * obj.units.k_length^2 / GapLength / (4 * pi * 1e-7);

            end

        end

        function torque = evalForceByArkkio(obj, mzName, GapLength)

            if obj.m.isTL3
                torque = obj.evalForceByArkkioTL3(mzName, GapLength);
            elseif obj.m.isTL6
                torque = obj.evalForceByArkkioTL6(mzName, GapLength);
            else
                error('Wrong element type');
            end

        end

        % force evaluation by Arrkio method using smoothed field data
        function torque = evalForceByArkkioSmoothTL3(obj, mzName, GapLength)

            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            bx = sum(obj.results.BxnSmooth(:,obj.m.ezi(:, mzptr.zi)),1)'/3;
            by = sum(obj.results.BynSmooth(:,obj.m.ezi(:, mzptr.zi)),1)'/3;
            b = [bx,by];
            p = mzptr.getCenterOfElements;
            b = ext_xy2rt(p, b);
            r = sqrt(sum(p.^2, 2));
            torque = mzptr.getAreaOfElements' * (b(:, 1) .* b(:, 2) .* r);
            torque = torque * obj.getDepth * obj.units.k_length^2 / GapLength / (4 * pi * 1e-7);

        end

        function torque = evalForceByArkkioSmooth(obj, mzName, GapLength)

            if obj.m.isTL3
                torque = obj.evalTorqueByArkkioSmoothTL3(mzName, GapLength);
            elseif obj.m.isTL6
                torque = obj.evalTorqueByArkkioTL6(mzName, GapLength);
            else
                error('Wrong element type');
            end

        end

        %% post-proccessing: magnitude plots
        function varargout = plotAmag(obj, varargin)

            % get figure, axes, and the number of contours
            [f,ax,Ncontour] = emdlab_flib_faxN(varargin{:});

            % plot Amag using patch function
            if obj.m.isTL3
                patch('faces', obj.m.cl(:, 1:3), 'vertices', obj.m.nodes, ...
                    'FaceVertexCData', obj.results.A, 'FaceColor', 'interp', ...
                    'EdgeColor', 'none', 'parent', ax);
            elseif obj.m.isTL6
                patch('faces', obj.m.cl(:, [1,4,2,5,3,6]), 'vertices', obj.m.nodes, ...
                    'FaceVertexCData', obj.results.A, 'FaceColor', 'interp', ...
                    'EdgeColor', 'none', 'parent', ax);
            else
                error('Wrong element type');
            end

            colormap(ax,jet(Ncontour));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Magnetic Vector Potential [Wb/m]';

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch('faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k', 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotFluxLines(obj, varargin)

            % get figure, axes, and the number of contours
            [f,ax,Ncontour] = emdlab_flib_faxN(varargin{:});

            % evaluation of contour lines
            if obj.m.isTL3
                cRange = linspace(min(obj.results.A), max(obj.results.A), Ncontour+2);
                c = tmzpc_contour_tl3(obj.m.cl, obj.m.nodes, obj.results.A, cRange(2:end-1));
                t = 1:size(c, 1);
                t = reshape(t, 2, [])';
            elseif obj.m.isTL6
                cRange = linspace(min(obj.results.A), max(obj.results.A), Ncontour+2);
                cl4 = [obj.m.cl(:,[1,4,6]);obj.m.cl(:,[4,2,5]);obj.m.cl(:,[6,4,5]);obj.m.cl(:,[6,5,3])];
                c = tmzpc_contour_tl3(cl4, obj.m.nodes, obj.results.A, cRange(2:end-1));
                t = 1:size(c, 1);
                t = reshape(t, 2, [])';
            else
                error('Wrong element type');
            end

            % plot mesh wireframe
            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch(ax, 'Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'k', 'EdgeColor', 'k', 'LineWidth', 0.6);

            % plot contour lines
            patch(ax, 'faces', t, 'vertices', c, 'edgecolor', 'b');

            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotMU(obj, varargin)

            if numel(varargin)
                [f,ax] = emdlab_flib_fax(varargin{1});
            else
                [f,ax] = emdlab_flib_fax();
            end

            % amplitude of the B at mesh points
            if numel(varargin)
                if isa(varargin{1}, 'matlab.graphics.axis.Axes')
                    ti = obj.m.getti(varargin{2:end});
                else
                    ti = obj.m.getti(varargin{:});
                end
            else
                ti = obj.m.getti;
            end
            ampB = (1./obj.edata.MagneticReluctivity(:,ti))/(4*pi*1e-7);
            ampB = repmat(ampB,3,1);
            %             ampB = sqrt(obj.results.Bxn(:,ti).^2 + obj.results.Byn(:,ti).^2);

            % plot using patch function
            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(ti, [1,2,3]))', ...
                'CData', ampB, 'FaceColor', 'interp', ...
                'EdgeColor', 'none', 'parent', ax);

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', [150,150,150]/255);
            colormap(ax,jet(15));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [tesla]';

            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotBmag(obj, varargin)

            if numel(varargin)
                [f,ax] = emdlab_flib_fax(varargin{1});
            else
                [f,ax] = emdlab_flib_fax();
            end

            % amplitude of the B at mesh points
            if numel(varargin)
                if isa(varargin{1}, 'matlab.graphics.axis.Axes')
                    ti = obj.m.getti(varargin{2:end});
                else
                    ti = obj.m.getti(varargin{:});
                end
            else
                ti = obj.m.getti;
            end
            ampB = sqrt(obj.results.Bxn(:,ti).^2 + obj.results.Byn(:,ti).^2);

            % plot using patch function
            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(ti, [1,2,3]))', ...
                'CData', ampB, 'FaceColor', 'interp', ...
                'EdgeColor', 'none', 'parent', ax);

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', [150,150,150]/255);
            colormap(ax,jet(15));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [tesla]';

            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotHmag(obj, varargin)

            % amplitude of the B at mesh points
            ti = obj.m.getti(varargin{:});
            ampH = sqrt(obj.results.Hxn(:,ti).^2 + obj.results.Hyn(:,ti).^2);

            % plot using patch function
            f = figure;
            ax = axes(f);
            f.Name = 'Magnetic flux density amplitude';

            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(ti, [1,2,3]))', ...
                'CData', ampH, 'FaceColor', 'interp', ...
                'EdgeColor', 'none', 'parent', ax);

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k');

            colormap(jet(15));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [tesla]';

            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = plotBmag3D(obj, varargin)

            % specefying zones
            ti = obj.m.getti(varargin{:});

            ampB = sqrt(sum(obj.results.Bxg(ti).^2 + obj.results.Byg(ti).^2, 2));
            % plot through patch
            f = figure;
            ax = axes(f);
            f.Name = '[Magnetic Flux Density Amplitude [Tesla]]';
            patch(ax, 'Faces', obj.m.cl(ti, 1:3), 'Vertices', [obj.m.nodes,ampB], ...
                'FaceVertexCData', ampB, 'FaceColor', 'flat', ...
                'EdgeColor', 'none');
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k');
            colormap(jet);
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [T]';
            %             cb.Limits(2) = 1.7;


            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = plotBmagF(obj, N, varargin)

            if nargin<2, N = 14; end

            obj.evalBn;

            % specefying zones
            ti = obj.m.getti(varargin{:});

            % calculate amplitude of B at mesh nodes
            ampB = sqrt(obj.results.BxnSmooth(:,ti).^2 + obj.results.BynSmooth(:,ti).^2);

            % plot using patch function
            f = figure;
            ax = axes(f);
            f.Name = 'Magnetic flux density with flux lines: smoothed';

            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(ti, [1,2,3]))', ...
                'CData', ampB, 'FaceColor', 'interp', 'EdgeColor', 'none', 'parent', ax);

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', [150,150,150]/255);
            colormap(jet(15));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [tesla]';

            % evaluation of contour lines
            if obj.m.isTL3
                cRange = linspace(min(obj.results.A), max(obj.results.A), N+2);
                c = tmzpc_contour_tl3(obj.m.cl, obj.m.nodes, obj.results.A, cRange(2:end-1));
                t = 1:size(c, 1);
                t = reshape(t, 2, [])';
            elseif obj.m.isTL6
                cRange = linspace(min(obj.results.A), max(obj.results.A), N+2);
                cl4 = [obj.m.cl(:,[1,4,6]);obj.m.cl(:,[4,2,5]);obj.m.cl(:,[6,4,5]);obj.m.cl(:,[6,5,3])];
                c = tmzpc_contour_tl3(cl4, obj.m.nodes, obj.results.A, cRange(2:end-1));
                t = 1:size(c, 1);
                t = reshape(t, 2, [])';
            else
                error('Wrong element type');
            end

            % plot contour lines
            patch(ax, 'faces', t, 'vertices', c, 'edgecolor', 'k');

            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotBmagFluxLines(obj, varargin)

            [f,ax,Ncontour] = emdlab_flib_faxN(varargin{:});

            obj.evalBn;

            % specefying zones
            ti = obj.m.getti();

            % calculate amplitude of B at mesh nodes
            ampB = sqrt(obj.results.BxnSmooth.^2 + obj.results.BynSmooth.^2);


            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(:, [1,2,3]))', ...
                'CData', ampB, 'FaceColor', 'interp', 'EdgeColor', 'none', 'parent', ax);

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', [150,150,150]/255);
            colormap(ax,jet(15));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [tesla]';

            % evaluation of contour lines
            if obj.m.isTL3
                cRange = linspace(min(obj.results.A), max(obj.results.A), Ncontour+2);
                c = tmzpc_contour_tl3(obj.m.cl, obj.m.nodes, obj.results.A, cRange(2:end-1));
                t = 1:size(c, 1);
                t = reshape(t, 2, [])';
            elseif obj.m.isTL6
                cRange = linspace(min(obj.results.A), max(obj.results.A), Ncontour+2);
                cl4 = [obj.m.cl(:,[1,4,6]);obj.m.cl(:,[4,2,5]);obj.m.cl(:,[6,4,5]);obj.m.cl(:,[6,5,3])];
                c = tmzpc_contour_tl3(cl4, obj.m.nodes, obj.results.A, cRange(2:end-1));
                t = 1:size(c, 1);
                t = reshape(t, 2, [])';
            else
                error('Wrong element type');
            end

            % plot contour lines
            patch(ax, 'faces', t, 'vertices', c, 'edgecolor', 'k');

            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = plotBmagFSkew(obj, ax, z, N, varargin)

            if nargin<4, N = 14; end

            obj.evalBn;

            % specefying zones
            ti = obj.m.getti(varargin{:});

            % calculate amplitude of B at mesh nodes
            ampB = sqrt(obj.results.BxnSmooth(:,ti).^2 + obj.results.BynSmooth(:,ti).^2);

            % plot using patch function
            f = ax.Parent;
            f.Name = 'Magnetic flux density with flux lines: smoothed skew';

            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            zdata = z * ones(obj.m.Nn,1);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(ti, [1,2,3]))', ...
                'ZData', zdata(obj.m.cl(ti, [1,2,3]))', ...
                'CData', ampB, 'FaceColor', 'interp', 'EdgeColor', 'none', 'parent', ax);

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', [obj.m.nodes,zdata], ...
                'EdgeColor', [150,150,150]/255);
            colormap(jet(15));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [tesla]';

            cRange = linspace(min(obj.results.A), max(obj.results.A), N+2);
            c = tmzpc_contour_tl3(obj.m.cl, obj.m.nodes, obj.results.A, cRange(2:end-1));
            t = 1:size(c, 1);
            t = reshape(t, 2, [])';
            patch(ax, 'faces', t, 'vertices', [c,z*ones(size(c,1),1)], 'edgecolor', 'k');

            %             axis off equal;
            %             zoom on;
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = plotBmagSmooth(obj, varargin)

            if numel(varargin)
                [f,ax] = emdlab_flib_fax(varargin{1});
            else
                [f,ax] = emdlab_flib_fax();
            end

            % amplitude of the B at mesh points
            if numel(varargin)
                if isa(varargin{1}, 'matlab.graphics.axis.Axes')
                    ti = obj.m.getti(varargin{2:end});
                else
                    ti = obj.m.getti(varargin{:});
                end
            else
                ti = obj.m.getti;
            end

            % calculate amplitude of B at mesh nodes
            ampB = sqrt(obj.results.BxnSmooth(:,ti).^2 + obj.results.BynSmooth(:,ti).^2);

            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(ti, [1,2,3]))', ...
                'CData', ampB, 'FaceColor', 'interp', 'EdgeColor', 'none', 'parent', ax);

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', [150,150,150]/255);
            colormap(jet(15));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [tesla]';

            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotHmagSmooth(obj, varargin)

            % calculate amplitude of B at mesh nodes
            obj.evalHn;
            ti = obj.m.getti(varargin{:});
            ampH = sqrt(obj.results.HxnSmooth(:,ti).^2 + obj.results.HynSmooth(:,ti).^2);

            % plot using patch function
            f = emdlab_r2d_mesh;
            ax = axes(f);
            f.Name = 'Magnetic field intensity amplitude: smoothed';

            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(ti, [1,2,3]))', ...
                'CData', ampH, 'FaceColor', 'interp', 'EdgeColor', 'none', 'parent', ax);

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k');

            colormap(jet(15));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Magnetic Field [A/m]';

            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        %% post-proccessing: vector plots
        function varargout = plotFQvecOnCenterOfElements(obj, FQname, varargin)

            if numel(varargin)
                [f,ax] = emdlab_flib_fax(varargin{1});
            else
                [f,ax] = emdlab_flib_fax();
            end

            ax.NextPlot = 'add';

            % get center of elements
            c = obj.m.getCenterOfElements;

            switch FQname
                case 'B'
                    obj.evalBn;
                    [FQx,FQy] = obj.getFQxFQyOnPoints(c(:,1),c(:,2),'B');

                case 'M'
                    FQx = obj.edata.MagnetizationX(1,:);
                    FQy = obj.edata.MagnetizationY(1,:);

                case 'H'
                    obj.evalHn;
                    [FQx,FQy] = obj.getFQxFQyOnPoints(c(:,1),c(:,2),'H');

                otherwise
                    error('Wrong field quantity.');
            end

            % mesh zone colors
            color = zeros(obj.m.Ne, 3);
            mzNames = fieldnames(obj.m.mzs);

            for i = 1:obj.m.Nmzs
                mzptr = obj.m.mzs.(mzNames{i});
                color(obj.m.ezi(:, mzptr.zi), :) = repmat(mzptr.color, mzptr.Ne, 1);
            end

            % plot mesh elements
            patch(ax, 'faces', obj.m.cl(:, 1:3), 'vertices', obj.m.nodes, ...
                'facecolor', 'flat', 'FaceVertexCData', color, 'EdgeColor', 'w', 'facealpha', 0.5);

            % plot vector field
            quiver(ax, c(:, 1), c(:, 2), FQx', FQy', 'color', 'k');

            % plot wireframe geometry
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', [130,130,130]/255);

            axis(ax, 'off', 'equal');
            zoom on;
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function plotBvecOnCenterOfElements(obj, varargin)
            obj.plotFQvecOnCenterOfElements('B', varargin{:});
        end

        function plotHvecOnCenterOfElements(obj, varargin)
            obj.plotFQvecOnCenterOfElements('H', varargin{:});
        end

        function plotMvecOnCenterOfElements(obj, varargin)
            obj.plotFQvecOnCenterOfElements('M', varargin{:});
        end

        function varargout = plotBvecOnNodes(obj)

            obj.evalBn;
            f = GraphicWindow;
            f.Name = 'Magnetic field on mesh nodes';
            h = guihandles(f);
            h.va.NextPlot = 'add';
            color = zeros(obj.m.Ne, 3);
            mzNames = fieldnames(obj.m.mzs);

            for i = 1:obj.m.Nmzs
                mzptr = obj.m.mzs.(mzNames{i});
                color(obj.m.ezi(:, mzptr.zi), :) = repmat(mzptr.color, mzptr.Ne, 1);
            end

            patch(h.va, 'faces', obj.m.cl(:, 1:3), 'vertices', obj.m.nodes, ...
                'facecolor', 'flat', 'FaceVertexCData', color, 'EdgeColor', 'w', 'facealpha', 1);
            quiver(h.va, obj.m.nodes(:, 1), obj.m.nodes(:, 2), obj.results.Bnx, obj.results.Bny, 'color', 'k');
            axis(h.va, 'off', 'equal');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function plotBvecOnGrid(obj, gridSize)

            if nargin < 2
                gridSize = 1;
            end

            tmp = min(obj.m.nodes);
            xmin = tmp(1);
            ymin = tmp(2);
            tmp = max(obj.m.nodes);
            xmax = tmp(1);
            ymax = tmp(2);
            [xg, yg] = meshgrid(linspace(xmin, xmax, ceil(xmax - xmin) / gridSize), ...
                linspace(ymin, ymax, ceil(ymax - ymin) / gridSize));
            xg = xg(:);
            yg = yg(:);
            e = pointLocation(triangulation(obj.m.cl, obj.m.nodes), [xg, yg]);
            index = ~ isnan(e);
            e = e(index);
            xg = xg(index);
            yg = yg(index);
            ah = setFigure();
            quiver(xg, yg, obj.results.Bxg(e), obj.results.Byg(e), 'color', 'b', 'parent', ah);
            axis off equal
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch('faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k', 'parent', ah);
            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end

        function plotBvecOnCircle(obj, varargin)

            [bx, by, x, y] = obj.getBxByOnCircle(varargin{:});
            close all
            quiver(x, y, bx, by);
            axis off equal;

        end

        function plotBrBtOnCircle(obj, varargin)

            [br, bt, t] = obj.getBrBtOnCircle(varargin{:});
            figure('Name', 'Br and Bt on circle')
            subplot(211)
            t = t * 180/pi;
            plot(t,br);
            set(gca,'Xlim',[0,360])
            xlabel('Mechanical Angle [deg]')
            ylabel('B_r [tesla]')
            title('Rdaial Flux Density Waveform');
            subplot(212)
            plot(t,bt)
            set(gca,'Xlim',[0,360])
            xlabel('Mechanical Angle [deg]')
            ylabel('B_t [tesla]')
            title('Tangentional Flux Density Waveform');
            zoom on;

        end

        function plotBrBtOnArc(obj, varargin)

            [br, bt, t] = obj.getBrBtOnArc(varargin{:});
            figure('Name', 'Br and Bt on circle')
            subplot(211)
            t = t * 180/pi;
            plot(t,br);
            set(gca,'Xlim',[varargin{4},varargin{5}]*180/pi)
            xlabel('Mechanical Angle [deg]')
            ylabel('B_r [tesla]')
            title('Rdaial Flux Density Waveform');
            subplot(212)
            plot(t,bt)
            set(gca,'Xlim',[varargin{4},varargin{5}]*180/pi)
            xlabel('Mechanical Angle [deg]')
            ylabel('B_t [tesla]')
            title('Tangentional Flux Density Waveform');
            zoom on;

        end

        function plotHrHtOnCircle(obj, varargin)

            [Hr, Ht, t] = obj.getHrHtOnCircle(varargin{:});
            figure('Name', 'Hr and Ht on circle')
            subplot(211)
            t = t * 180/pi;
            plot(t,Hr/1e3);
            set(gca,'Xlim',[0,360])
            xlabel('Mechanical Angle [deg]')
            ylabel('H_r [kA/m]')
            title('Rdaial Magnetic Field Waveform');
            subplot(212)
            plot(t,Ht/1e3)
            set(gca,'Xlim',[0,360])
            xlabel('Mechanical Angle [deg]')
            ylabel('H_t [kA/m]')
            title('Tangentional Magnetic Field Waveform');
            zoom on;

        end

        function plotBrFFT(obj, poles, Nharmonics, varargin)

            obj.evalBn;
            index = obj.m.getnIndexOnCircle(varargin{:});
            p = obj.m.nodes(index, :);
            b = [obj.Bn.x(index), obj.Bn.y(index)];
            pAngle = atan_02pi(p);
            [pAngle, index] = sort(pAngle);
            p = p(index, :);
            b = b(index, :);
            b = ext_xy2rt(p, b);
            figure('Name', '[@EMDLab] Middle Air Gap Br and Bt', 'NumberTitle', 'off', ...
                'WindowStyle', 'modal')
            subplot(211)
            plot(pAngle * 180 / pi, b(:, 1))
            set(gca, 'Xlim', [pAngle(1) * 180 / pi, pAngle(end) * 180 / pi])
            xlabel('Mechanical Angle [Degree]')
            ylabel('B_r [Tesla]')
            title(['Rdaial Flux Density, Mean = ', num2str(mean(abs(b(:, 1))))])
            hold all;
            an = (1 / pi) * int_trap(full(repmat(b(:, 1)', Nharmonics, 1) .* cos(poles * sparse(1:Nharmonics, 1:Nharmonics, 1:Nharmonics) * repmat(pAngle', Nharmonics, 1))), pAngle');
            bn = (1 / pi) * int_trap(full(repmat(b(:, 1)', Nharmonics, 1) .* sin(poles * sparse(1:Nharmonics, 1:Nharmonics, 1:Nharmonics) * repmat(pAngle', Nharmonics, 1))), pAngle');
            hn = sqrt(an.^2 + bn.^2);
            y = zeros(length(pAngle), 1);

            for i = 1:Nharmonics
                y = y + an(i) * cos(poles * i * pAngle) + bn(i) * sin(poles * i * pAngle);
                plot(pAngle * 180 / pi, y, '-.');
            end

            subplot(212)
            stem(1:Nharmonics, hn)
            set(gca, 'Xlim', [0, Nharmonics + 1])
            xlabel('Harmonic Index')
            ylabel('Harmonic Amplitude')
            title(['Harmonic Decomposition of Radial Flux, H[1] = ', num2str(hn(1))]);

        end

        %% post-proccessing: interpolation functions
        % get triangle index (ti) of points specified by theri (x,y) coordinates
        function ti = getPointsLocation(obj, x, y)

            if (~isvector(x)) || (~isvector(y))
                MException('', 'x and y inputs must be row vectors.');
            end

            if length(x) ~= length(y)
                MException('', 'Length of x must be the same as length of y.');
            end

            if isrow(x), x = x'; end
            if isrow(y), y = y'; end

            % find index of triangle containing the (x,y) point
            if obj.m.isTL3
                ti = pointLocation(triangulation(obj.m.cl, obj.m.nodes), [x,y]);
            elseif obj.m.isTL6
                ti = pointLocation(triangulation(obj.m.cl(:,1:3), obj.m.nodes(1:end-size(obj.m.edges,1),:)), [x,y]);
            else
                error('Wrong element type');
            end

        end

        % get the value of magnetic vector potential Az on points
        function Az = getAzOnPoints(obj, x, y)

            % find triangle indicies related to point locations
            ti = obj.getPointsLocation(x, y);

            % care points out of the mesh
            index = isnan(ti);
            ti(index) = 1;

            % interpolate the value of Az on points
            if obj.m.isTL3
                Az = emdlab_m2d_tl3_interpA(obj.m.cl, obj.m.nodes, obj.results.A(obj.m.cl'), obj.m.JIT, ti, x, y);
            elseif obj.m.isTL6
                Az = emdlab_m2d_tl6_interpA(obj.m.cl, obj.m.nodes, obj.results.A(obj.m.cl'), obj.m.JIT, ti, x, y);
            else
                error('Wrong element type');
            end

            Az(index) = NaN;

            if iscolumn(Az), Az = Az'; end

        end

        % get x- and y-components of a field quantity: B, BSmooth, H, HSmooth
        function [FQx, FQy] = getFQxFQyOnPointsTL3(obj, x, y, fieldQuantity)

            % find index of triangle containing the (x,y) point
            ti = obj.getPointsLocation(x, y);

            % care points out of the mesh
            index = isnan(ti);
            ti(index) = 1;

            switch fieldQuantity

                case 'B'
                    FQx = obj.results.Bxg(ti);
                    FQy = obj.results.Byg(ti);

                case 'BSmooth'
                    [FQx, FQy] = emdlab_m2d_tl3_interpFQSmooth(obj.m.cl, obj.m.nodes, obj.results.Bxn, obj.results.Byn, obj.m.JIT, ti, x, y);

                case 'H'
                    FQx = obj.results.Hxg(ti);
                    FQy = obj.results.Hyg(ti);

                case 'HSmooth'
                    [FQx, FQy] = emdlab_m2d_tl3_interpFQSmooth(obj.m.cl, obj.m.nodes, obj.results.Hxn, obj.results.Hyn, obj.m.JIT, ti, x, y);

                otherwise
                    error('Wrong field quantity.');
            end

            FQx(index) = NaN;
            FQy(index) = NaN;

            if iscolumn(FQx), FQx = FQx'; end
            if iscolumn(FQy), FQy = FQy'; end

        end

        function [FQx, FQy] = getFQxFQyOnPointsTL6(obj, x, y, fieldQuantity)

            % find index of triangle containing the (x,y) point
            ti = obj.getPointsLocation(x, y);

            % care points out of the mesh
            index = isnan(ti);
            ti(index) = 1;

            switch fieldQuantity

                case 'B'
                    [FQx, FQy] = emdlab_m2d_tl6_interpB(obj.m.cl, obj.m.nodes, obj.results.A, obj.m.JIT, ti, x, y);

                case 'BSmooth'
                    [FQx, FQy] = emdlab_m2d_tl6_interpBSmooth(obj.m.cl, obj.m.nodes, obj.results.BxnSmooth, obj.results.BynSmooth, obj.m.JIT, ti, x, y);

                case 'H'
                    [FQx, FQy] = emdlab_m2d_tl6_interpBSmooth(obj.m.cl, obj.m.nodes, obj.results.Hxn, obj.results.Hyn, obj.m.JIT, ti, x, y);

                case 'HSmooth'
                    [FQx, FQy] = emdlab_m2d_tl6_interpBSmooth(obj.m.cl, obj.m.nodes, obj.results.Hxn, obj.results.Hyn, obj.m.JIT, ti, x, y);

                otherwise
                    error('Wrong field quantity.');
            end

            FQx(index) = NaN;
            FQy(index) = NaN;

            if iscolumn(FQx), FQx = FQx'; end
            if iscolumn(FQy), FQy = FQy'; end

        end

        function [FQx, FQy] = getFQxFQyOnPoints(obj, x, y, fieldQuantity)
            if obj.m.isTL3
                [FQx, FQy] = obj.getFQxFQyOnPointsTL3(x, y, fieldQuantity);
            elseif obj.m.isTL6
                [FQx, FQy] = obj.getFQxFQyOnPointsTL6(x, y, fieldQuantity);
            else
                error('Wrong element type');
            end
        end

        % get the value of Bx and By on specified points
        function [Bx, By] = getBxByOnPoints(obj, x, y, smoothFlag)

            % set default value of smoothFlag
            if nargin < 4, smoothFlag = true; end

            if smoothFlag
                [Bx, By] = getFQxFQyOnPoints(obj, x, y, 'BSmooth');
            else
                [Bx, By] = getFQxFQyOnPoints(obj, x, y, 'B');
            end

        end

        % get the value of Hx and Hy on specified points
        function [Hx, Hy] = getHxHyOnPoints(obj, x, y, smoothFlag)

            % set default value of smoothFlag
            if nargin < 4, smoothFlag = true; end

            if smoothFlag
                [Hx, Hy] = getFQxByFQnPoints(obj, x, y, 'HSmooth');
            else
                [Hx, Hy] = getFQxByFQnPoints(obj, x, y, 'H');
            end

        end

        % get x- and y-components of a field quantity on a segment
        function [FQx, FQy, x, y] = getFQxFQyOnSegment(obj, x1, y1, x2, y2, N, fieldQuantity)

            % sample points
            x = linspace(x1, x2, N);
            y = linspace(y1, y2, N);
            [FQx, FQy] = obj.getFQxFQyOnPoints(x, y, fieldQuantity);

        end

        % get the value of Bx and By on points that are on a segment
        function [Bx, By, x, y] = getBxByOnSegment(obj, x1, y1, x2, y2, N)

            % set default value of N
            if nargin < 6, N = 1000; end
            [Bx, By, x, y] = obj.getFQxFQyOnSegment(x1, y1, x2, y2, N, 'BSmooth');

        end

        % get the value of Hx and Hy on points that are on a segment
        function [Hx, Hy, x, y] = getHxHyOnSegment(obj, x1, y1, x2, y2, N)

            % set default value of N
            if nargin < 6, N = 1000; end
            [Hx, Hy, x, y] = obj.getFQxFQyOnSegment(x1, y1, x2, y2, N, 'HSmooth');

        end

        % get x- and y-components of a field quantity on a circle
        function [FQx, FQy, x, y] = getFQxFQyOnCircle(obj, xc, yc, r, N, fieldQuantity)

            % sample points
            t = linspace(0,2*pi,N);
            x = xc + r * cos(t);
            y = yc + r * sin(t);
            [FQx, FQy] = obj.getFQxFQyOnPoints(x, y, fieldQuantity);

        end

        % get x- and y-components of a field quantity on an arc
        function [FQx, FQy, x, y] = getFQxFQyOnArc(obj, xc, yc, r, theta1, theta2, N, fieldQuantity)

            % sample points
            t = linspace(theta1, theta2, N);
            x = xc + r * cos(t);
            y = yc + r * sin(t);
            [FQx, FQy] = obj.getFQxFQyOnPoints(x, y, fieldQuantity);

        end

        % get the value of Bx and By on points that are on a circle
        function [Bx, By, x, y] = getBxByOnCircle(obj, xc, yc, r, N)

            % set default value of N
            if nargin < 5, N = 1000; end
            [Bx, By, x, y] = obj.getFQxFQyOnCircle(xc, yc, r, N, 'BSmooth');

        end

        % get the value of Hx and Hy on points that are on a circle
        function [Hx, Hy, x, y] = getHxHyOnCircle(obj, xc, yc, r, N)

            % set default value of N
            if nargin < 5, N = 1000; end
            [Hx, Hy, x, y] = obj.getFQxFQyOnCircle(xc, yc, r, N, 'HSmooth');

        end

        % get the value of Bx and By on points that are on a circle
        function [bx, by, x, y] = getBxByOnArc(obj, xc, yc, r, a1, a2, N)

            % set default value of N
            if nargin < 7, N = 1000; end

            % sample points
            t = linspace(a1,a2,N);
            x = xc + r * cos(t);
            y = yc + r * sin(t);
            [bx, by] = obj.getBxByOnPoints(x,y);

        end

        % get the value of FQr and FQt on points that are on a circle
        function [FQr, FQt, t] = getFQrFQtOnCircle(obj, xc, yc, r, N, fieldQuantity)

            [FQx, FQy, x, y] = obj.getFQxFQyOnCircle(xc, yc, r, N, fieldQuantity);
            t = linspace(0,2*pi,N);

            % unit vector in radial direction
            uRx = x - xc;
            uRy = y - yc;
            uRMagnitude = sqrt(uRx.^2 + uRy.^2);
            uRx = uRx./uRMagnitude;
            uRy = uRy./uRMagnitude;

            % perform inner products to calculate br and bt
            FQr = FQx.*uRx + FQy.*uRy;
            FQt = -FQx.*uRy + FQy.*uRx;

        end

        % get the value of FQr and FQt on points that are on an arc
        function [FQr, FQt, t] = getFQrFQtOnArc(obj, xc, yc, r, theta1, theta2, N, fieldQuantity)

            [FQx, FQy, x, y] = obj.getFQxFQyOnArc(xc, yc, r, theta1, theta2, N, fieldQuantity);
            t = linspace(theta1,theta2,N);

            % unit vector in radial direction
            uRx = x - xc;
            uRy = y - yc;
            uRMagnitude = sqrt(uRx.^2 + uRy.^2);
            uRx = uRx./uRMagnitude;
            uRy = uRy./uRMagnitude;

            % perform inner products to calculate br and bt
            FQr = FQx.*uRx + FQy.*uRy;
            FQt = -FQx.*uRy + FQy.*uRx;

        end

        % get the value of Br and Bt on points that are on a circle
        function [Br, Bt, t] = getBrBtOnCircle(obj, xc, yc, r, N)

            % set default value of N
            if nargin < 5, N = 1000; end
            [Br, Bt, t] = obj.getFQrFQtOnCircle(xc, yc, r, N, 'BSmooth');

        end

        % get the value of Br and Bt on points that are on a circle
        function [Br, Bt, t] = getBrBtOnArc(obj, xc, yc, r, theta1, theta2, N)

            % set default value of N
            if nargin < 7, N = 1000; end
            [Br, Bt, t] = obj.getFQrFQtOnArc(xc, yc, r, theta1, theta2, N, 'BSmooth');

        end

        % get the value of Br and Bt on points that are on a circle
        function [Hr, Ht, t] = getHrHtOnCircle(obj, xc, yc, r, N)

            % set default value of N
            if nargin < 5, N = 1000; end
            [Hr, Ht, t] = obj.getFQrFQtOnCircle(xc, yc, r, N, 'HSmooth');

        end

    end

    % solver methods
    methods

        function y = getTotalEnergy(obj)
            y = obj.evalTotalEnergyCoenergy;
        end

        function y = getTotalCoenergy(obj)
            [~,y] = obj.evalTotalEnergyCoenergy;
        end

        function y = getDepth(obj)
            y = obj.depth * obj.units.k_length;
        end

        function y = getMeshZoneAc(obj, mzName)
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            y = obj.results.A(mzptr.l2g);
            y = y(mzptr.cl);
            y = mean(y, 2);
        end

        function f = gui(obj)

            f = figure('Position',[0,0,850,600], 'Visible','off', 'Name', 'https://github.com/EMDLAB-Package/emdlab-win64', 'DockControls', 'off');
            movegui(f,'center');
            tgp = uitabgroup(f, 'Tag', 'MainTabGroup');

            t = uitab(tgp,'Title','Geometry');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.m.showg(ax);

            t = uitab(tgp,'Title','Mesh');

            tgp_sub = uitabgroup(t);

            t = uitab(tgp_sub,'Title','Mesh Zones');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.m.showmzs(ax);

            t = uitab(tgp_sub,'Title','Global Mesh');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.m.showm(ax);

            t = uitab(tgp_sub,'Title','Wireframe Mesh');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.m.showwf(ax);

            t = uitab(tgp_sub,'Title','Free Boundary');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.m.showfb(ax);

            t = uitab(tgp_sub,'Title','Mesh Information');
            tbl = uitable(t, 'units', 'normalized', 'position', [0,0,1,1], 'fontSize',12);

            mzNames = obj.m.getMeshZoneNames;

            tbl.ColumnName = {'   Mesh Zone Name   ','   Number of Elements   '};
            tbl.Data{1,1} = 'All Mesh Zones';
            tbl.Data{1,2} = obj.m.Ne;
            for i = 1:obj.m.Nmzs
                tbl.Data{i+1,1} = char(mzNames(i));
                tbl.Data{i+1,2} = obj.m.mzs.(mzNames(i)).Ne;
            end
            tbl.ColumnWidth = {'auto'};

            t = uitab(tgp,'Title','Coils');
            tgp_sub = uitabgroup(t);

            t = uitab(tgp_sub,'Title','#All coil arms');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.showAllCoilArms(ax);

            coilsNames = fieldnames(obj.coils);
            for i = 1:obj.Ncoils
                t = uitab(tgp_sub,'Title',coilsNames{i});
                ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
                obj.showCoil(coilsNames{i}, ax);
            end

            t = uitab(tgp_sub,'Title','Flux Linkages');
            tbl = uitable(t, 'units', 'normalized', 'position', [0,0,1,1], 'fontSize',12);

            tbl.ColumnName = {'   Coil Name   ','   Flux Linkage [Wb]   '};
            for i = 1:obj.Ncoils
                tbl.Data{i,1} = coilsNames{i};
                tbl.Data{i,2} = obj.evalCoilFluxLinkage(coilsNames{i});
            end
            tbl.ColumnWidth = {'auto'};

            t = uitab(tgp,'Title','Boundary Conditions');
            tgp_sub = uitabgroup(t);

            bcNames = fieldnames(obj.bcs.dirichlet);
            for i = 1:numel(bcNames)
                t = uitab(tgp_sub,'Title','Vector Potential #'+string(i));
                ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
                obj.m.showNodes(ax,obj.bcs.dirichlet.(bcNames{i}).index);
            end

            bcNames = fieldnames(obj.bcs.evenPeriodic);
            for i = 1:numel(bcNames)
                t = uitab(tgp_sub,'Title','Even Periodic #'+string(i));
                ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
                obj.m.showNodes(ax,obj.bcs.evenPeriodic.(bcNames{i}).mIndex, obj.bcs.evenPeriodic.(bcNames{i}).sIndex);
            end

            bcNames = fieldnames(obj.bcs.oddPeriodic);
            for i = 1:numel(bcNames)
                t = uitab(tgp_sub,'Title','Odd Periodic #'+string(i));
                ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
                obj.m.showNodes(ax,obj.bcs.oddPeriodic.(bcNames{i}).mIndex, obj.bcs.oddPeriodic.(bcNames{i}).sIndex);
            end

            t = uitab(tgp,'Title','Field Plots');

            tgp_sub = uitabgroup(t);

            t = uitab(tgp_sub,'Title','Flux Density');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.plotBmag(ax);

            t = uitab(tgp_sub,'Title','Flux Density Smoothed');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.plotBmagSmooth(ax);

            t = uitab(tgp_sub,'Title','Flux Lines');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.plotFluxLines(ax);

            t = uitab(tgp_sub,'Title','Flux Density & Flux Lines');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.plotBmagFluxLines(ax);
            colormap(ax,"parula");

            t = uitab(tgp_sub,'Title','Flux Density Field');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.plotBvecOnCenterOfElements(ax);

            t = uitab(tgp_sub,'Title','Magnetization');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.plotMvecOnCenterOfElements(ax);

            t = uitab(tgp_sub,'Title','Vector Potential');
            ax = axes(t, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            obj.plotAmag(ax);

        end

    end

    methods (Access = protected)

        function makeFalse_isElementDataAssigned(obj)
            obj.isElementDataAssigned = false;
            obj.isBnEvaluated = false;
        end

    end

end

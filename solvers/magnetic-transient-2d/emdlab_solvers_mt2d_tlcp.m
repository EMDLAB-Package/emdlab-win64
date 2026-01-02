% EMDLAB: Electrical Machines Design Laboratory
% magnetic-transient two-dimensional tlcp
% triangular lagrangian elements common properties

classdef emdlab_solvers_mt2d_tlcp < handle & matlab.mixin.Copyable

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

        % cages
        cages (1,1) struct;

        % start connections
        starConnections (1,1) struct;

        % an structure to hold matrices
        mtcs (1,1) struct;

        % moving regions
        movingRegions (1,1) struct = struct();

        % solver index: is used for parallel computing
        solverIndex = 1;

    end

    properties (SetAccess = protected)

        % simulation time
        simTime (1,:) double = 0;

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
        isSolvedForInitialConditions (1,1) logical = false;

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

        % number of cages
        Ncages (1,1) double;

        % number of cages
        NstarConnections (1,1) double;

        % number of motions
        NmovingRegions (1,1) double;

    end

    methods
        %% constructor and destructor
        function obj = emdlab_solvers_mt2d_tlcp()

            % default values
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

        function y = get.Ncages(obj)
            y = numel(fieldnames(obj.cages));
        end

        function y = get.NstarConnections(obj)
            y = numel(fieldnames(obj.starConnections));
        end

        function y = get.NmovingRegions(obj)
            y = numel(fieldnames(obj.movingRegions));
        end

        %% solver properties for mesh zones
        % set default properties of a mesh zone
        function setdp(obj, mzName)

            obj.m.mzs.(mzName).props.isCoilArm = false;
            obj.m.mzs.(mzName).props.isMagnetized = false;
            obj.m.mzs.(mzName).props.isEddyZone = false;
            obj.m.mzs.(mzName).props.isMoving = false;
            obj.m.mzs.(mzName).props.isCoreLossActivated = false;
            obj.makeFalse_isElementDataAssigned;

        end

        function addmz(obj, mzName, mzValue)
            mzName = rmspaces(mzName);
            obj.m.addmz(mzName, mzValue);
            obj.setdp(mzName);
        end

        function obj = removemz(obj, mzName)
            obj.m.removemz(mzName);
        end

        function rotateMeshZone(obj, mzName, varargin)
            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.m.rotateMeshZone(mzName, varargin{:});

            if obj.m.mzs.(mzName).props.isMagnetized
                obj.m.mzs.(mzName).props.magnetization.rotate(varargin{:});
            end

        end

        function rotateMeshZones(obj, meshZoneNames, varargin)

            tmp = zeros([],1);
            for i = 1:numel(meshZoneNames)
                mzptr = obj.m.mzs.(meshZoneNames(i));

                if mzptr.props.isMoving
                    obj.rotateMeshZone(meshZoneNames(i), varargin{:});
                    tmp = [tmp;mzptr.l2g];
                    [mzptr.props.JIT,~] = emdlab_m2d_tl3_evalJIT(mzptr.cl, mzptr.nodes);
                end

            end
            tmp = unique(tmp);
            obj.m.nodes(tmp,:) = emdlab_g2d_rotatePoints(obj.m.nodes(tmp,:),varargin{:});           

        end

        %% excitation definitions => electrical or mechanical
        % define a new coil
        function coilName = checkCoilExistence(obj, coilName)

            if ~isfield(obj.coils, coilName)
                throw(MException('', ['Coil with name [', coilName, '] does not exist.']));
            end

        end

        function coilName = checkCoilNonExistence(obj, coilName)

            if isfield(obj.coils, coilName)
                throw(MException('', ['Another coil with name [', coilName, '] already exist.']));
            end

        end

        function defineCoil(obj, coilName, fedType, eddyType)

            if nargin < 3, fedType = 'voltage'; end
            if nargin < 4, eddyType = 'stranded'; end
            fedType = lower(erase(fedType, ' '));
            eddyType = lower(erase(eddyType, ' '));
            coilName = obj.checkCoilNonExistence(coilName);
            obj.coils.(coilName) = emdlab_solvers_mt2d_coil();
            obj.coils.(coilName).fedType = fedType;
            obj.coils.(coilName).eddyType = eddyType;
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

                if strcmpi(cptr.eddyType, 'solid')
                    obj.m.mzs.(mzName(i)).props.isEddyZone = true;
                end

            end

        end

        % this function add a mesh zone to a coil and make it as a coil arm
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

                if strcmpi(cptr.eddyType, 'solid')
                    obj.m.mzs.(mzName(i)).props.isEddyZone = true;
                end

            end

        end

        % this function add a mesh zone to a coil and make it as a coil arm
%         function addMeshZone2Coil(obj, coilName, mzName, turns, direction, kfill)
% 
%             % default arguments
%             if nargin < 4, turns = 1; end
%             if nargin < 5, direction = 1; end
%             if nargin < 6, kfill = 1; end
% 
%             coilName = obj.checkCoilExistence(coilName);
%             mzName = obj.m.checkMeshZoneExistence(mzName);
% 
%             if obj.m.mzs.(mzName).props.isCoilArm
%                 error('Specified mesh zone is already defined as a coil arm.');
%             end
% 
%             % get coil pointer
%             cptr = obj.coils.(coilName);
% 
%             if strcmpi(cptr.eddyType, 'solid')
%                 if turns ~= 1
%                     error('The number of turns for solid coil arm must be one.');
%                 end
%                 if kfill ~= 1
%                     error('The coil arm fill factor for solid coil arm must be one.');
%                 end
%             end
% 
%             if ~ismember(direction, [-1,1])
%                 error('The coil arm reference direction must be <1> or <-1>.');
%             end
% 
%             if kfill > 1
%                 error('The coil fill factor must be lower than or equal to one.');
%             end
% 
%             cptr.addCoilArm(mzName, direction);
%             obj.m.mzs.(mzName).props.turns = turns;
%             obj.m.mzs.(mzName).props.direction = direction;
%             obj.m.mzs.(mzName).props.kfill = kfill;
%             obj.m.mzs.(mzName).props.isCoilArm = true;
%             obj.coilArms(end+1) = mzName;
%             obj.m.mzs.(mzName).props.cai = obj.NcoilArms;
% 
%             if strcmpi(cptr.eddyType, 'solid')
%                 obj.m.mzs.(mzName).props.isEddyZone = true;
%             end
% 
%         end

        % check if all defined coils have coil arms
        function checkCoils(obj)

            % coils with zeros number of coil arms
            coilNames = fieldnames(obj.coils);
            for i = 1:obj.Ncoils

                if obj.coils.(coilNames{i}).NcoilArms == 0
                    error(['Coil <', coilNames{i}, '> does not have any coil arm.']);
                end

            end

            % calculate DC resistance of the coils
            for i = 1:obj.Ncoils

                % get coil pointer
                cptr = obj.coils.(coilNames{i});
                cptr.Rdc = 0;

                for j = 1:cptr.NcoilArms

                    % pointer to coil arm
                    mzptr = obj.m.mzs.(cptr.coilArms(j));
                    
                    % calculate DC resistance of coil arms
                    mzptr.props.Rdc = obj.getDepth * mzptr.props.turns^2 /(obj.m.mts.(mzptr.material).ElectricConductivity.value * ...
                        mzptr.getArea * mzptr.props.kfill * obj.units.k_length^2);
                    cptr.Rdc = cptr.Rdc + mzptr.props.Rdc;
                    
                end
 
            end

        end

        function computeCoilsRdc (obj)
            obj.checkCoils;
             % coils with zeros number of coil arms
            coilNames = fieldnames(obj.coils);
            for i = 1:obj.Ncoils
                disp("#" + string(i) + ": " + coilNames{i} + " -> Rdc = " + num2str(obj.coils.(coilNames{i}).Rdc));
            end
        end

        function setCoilVoltage(obj, coilName, value)
            coilName = obj.checkCoilExistence(coilName);
            obj.coils.(coilName).setVoltage(value);
        end

        function setCoilCurrent(obj, coilName, value)
            coilName = obj.checkCoilExistence(coilName);
            obj.coils.(coilName).setCurrent(value);
        end

        function setCoilInitialCurrent(obj, coilName, value)
            coilName = obj.checkCoilExistence(coilName);
            obj.coils.(coilName).setInitialCurrent(value);
        end

        function forcePeriodicCoilVoltages(obj)
            for cName = string(fieldnames(obj.coils)')
                cptr = obj.coils.(cName);
                cptr.inducedVoltage(1) = cptr.inducedVoltage(end);
                cptr.voltage(1) = cptr.voltage(end);
            end
        end

        % cage definition
        function cageName = checkCageExistence(obj, cageName)

            if ~isfield(obj.cages, cageName)
                throw(MException('', ['Cage with name [', cageName, '] does not exist.']));
            end

        end

        function cageName = checkCageNonExistence(obj, cageName)

            if isfield(obj.cages, cageName)
                throw(MException('', ['Another cage with name [', cageName, '] already exist.']));
            end

        end

        function defineCage(obj, cageName, meshZoneNames, Re, Le, cageType)

            if ~isstring(meshZoneNames) || (numel(meshZoneNames)<3)
                error('meshZoneNames must be a string vector, the number of mesh zones must be higher than 2.');
            end

            for i = 1:numel(meshZoneNames)
                obj.m.checkMeshZoneExistence(meshZoneNames(i));
            end

            if nargin < 4, Re = 0; end
            if nargin < 5, Le = 0; end
            if nargin < 6, cageType = 'even'; end
            cageType = lower(erase(cageType, ' '));
            cageName = obj.checkCageNonExistence(cageName);
            obj.cages.(cageName) = emdlab_solvers_mt2d_cage();

            for i = 1:numel(meshZoneNames)
                obj.defineCoil(meshZoneNames(i), 'voltage', 'solid');
                obj.addMeshZone2Coil(meshZoneNames(i), meshZoneNames(i), 1, 1, 1);
                obj.coils.(meshZoneNames(i)).isCageMember = true;

                % store important indices
                if i == 1
                    obj.cages.(cageName).caiStart = obj.NcoilArms;
                    obj.cages.(cageName).ciStart = obj.Ncoils;
                elseif i == numel(meshZoneNames)
                    obj.cages.(cageName).caiEnd = obj.NcoilArms;
                    obj.cages.(cageName).ciEnd = obj.Ncoils;
                end
            end

            obj.cages.(cageName).coilArms = meshZoneNames;
            obj.cages.(cageName).Re = Re;
            obj.cages.(cageName).Le = Le;
            obj.cages.(cageName).type = cageType;
            obj.cages.(cageName).type;

        end

         % cage definition
         function starConnectionName = checkStarConnectionExistence(obj, starConnectionName)

            if ~isfield(obj.starConnections, starConnectionName)
                throw(MException('', ['Star connection with name [', starConnectionName, '] does not exist.']));
            end

        end

        function starConnectionName = checkStarConnectionNonExistence(obj, starConnectionName)

            if isfield(obj.starConnections, starConnectionName)
                throw(MException('', ['Another Star connection with name [', starConnectionName, '] already exist.']));
            end

        end

        function defineStarConnection(obj, starConnectionName, varargin)

            starConnectionName = obj.checkStarConnectionNonExistence(starConnectionName);
            obj.starConnections.(starConnectionName) = emdlab_solvers_mt2d_star();
            obj.starConnections.(starConnectionName).sci = obj.NstarConnections;

            if all(cellfun(@ischar, varargin))
                coilNames = string(varargin);
            else
                coilNames = string([varargin{:}]);
            end
            
            % loop over coil names
            cIndicies = zeros(1,numel(coilNames));
            for i = 1:numel(coilNames)
                coilName = obj.checkCoilExistence(coilNames(i));
                cptr = obj.coils.(coilName);
                if cptr.isCurrentFed || cptr.isSolid
                    error('All star connection coils must be stranded voltage fed coils.');
                end
                cptr.isStarConnectionMember = true;
                cIndicies(i) = cptr.ci;
            end

            obj.starConnections.(starConnectionName).coilArms = coilNames;
            obj.starConnections.(starConnectionName).ci = cIndicies;

        end

        function setMagnetization(obj, mzName, varargin)
            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.m.mzs.(mzName).props.isMagnetized = true;

            if nargin == 3 && isa(varargin{1}, 'mdMagnetization')
                obj.m.mzs.(mzName).props.magnetization = varargin{1};
            else
                obj.m.mzs.(mzName).props.magnetization = emdlab_solvers_mt2d_magnetization(varargin{:});
            end

            % change states
            obj.makeFalse_isElementDataAssigned;
        end

        % activate eddy effect calculation for specfied mesh zones
        function activateEddyEffect(obj, varargin)
        
            % loop over inputs
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    mzName = obj.m.checkMeshZoneExistence(varargin{i});
                    obj.m.mzs.(mzName).props.isEddyZone = true;
                elseif isvector(varargin{i}) && isstring(varargin{i})
                    for j = 1:numel(varargin{i})
                        mzName = char(varargin{i}(j));
                        mzName = obj.m.checkMeshZoneExistence(mzName);
                        obj.m.mzs.(mzName).props.isEddyZone = true;
                    end
                else
                    error('The input class type must be <char> or <string>.');
                end
            end

        end

        % activate eddy effect calculation for specfied mesh zones
        function deactivateEddyEffect(obj, varargin)

            % loop over inputs
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    mzName = obj.m.checkMeshZoneExistence(varargin{i});
                    obj.m.mzs.(mzName).props.isEddyZone = false;
                elseif isvector(varargin{i}) && isstring(varargin{i})
                    for j = 1:numel(varargin{i})
                        mzName = char(varargin{i}(j));
                        mzName = obj.m.checkMeshZoneExistence(mzName);
                        obj.m.mzs.(mzName).props.isEddyZone = false;
                    end
                else
                    error('The input class type must be <char> or <string>.');
                end
            end

        end

        function [ySum,y] = evalSolidLoss(obj)
            [ySum,y] = emdlab_flib_calculateMECLM1(obj);
        end

        % cage definition
        function movingRegionName = checkMovingRegionExistence(obj, movingRegionName)

            if ~isfield(obj.movingRegions, movingRegionName)
                throw(MException('', ['Moving region with name [', movingRegionName, '] does not exist.']));
            end

        end

        function movingRegionName = checkMovingRegionNonExistence(obj, movingRegionName)

            if isfield(obj.movingRegions, movingRegionName)
                throw(MException('', ['Another moving region with name [', movingRegionName, '] already exist.']));
            end

        end

        function defineMovingRegion(obj, movingRegionName, meshZoneNames, interfaceMeshZone)

            for i = 1:numel(meshZoneNames)
                obj.m.checkMeshZoneExistence(meshZoneNames(i));
            end

            obj.m.checkMeshZoneExistence(interfaceMeshZone);
            if ~obj.m.mzs.(interfaceMeshZone).props.isInterface
                error("Mesh zone" + interfaceMeshZone + "is not defined as an interface mesh zone.");
            end

            movingRegionName = obj.checkMovingRegionNonExistence(movingRegionName);
            obj.movingRegions.(movingRegionName) = emdlab_solvers_mt2d_motion();
            obj.movingRegions.(movingRegionName).meshZones = meshZoneNames;
            obj.movingRegions.(movingRegionName).mi = obj.NmovingRegions;
            obj.movingRegions.(movingRegionName).interface = interfaceMeshZone;
            obj.setMoving(movingRegionName, meshZoneNames);

        end

        function setMoving(obj, movingRegionName, varargin)
            % loop over inputs
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    mzName = obj.m.checkMeshZoneExistence(varargin{i});
                    obj.m.mzs.(mzName).props.isMoving = true;
                    obj.m.mzs.(mzName).props.movingRegionName = movingRegionName;
                elseif isvector(varargin{i}) && isstring(varargin{i})
                    for j = 1:numel(varargin{i})
                        mzName = char(varargin{i}(j));
                        mzName = obj.m.checkMeshZoneExistence(mzName);
                        obj.m.mzs.(mzName).props.isMoving = true;
                        obj.m.mzs.(mzName).props.movingRegionName = movingRegionName;
                    end
                else
                    error('The input class type must be <char> or <string>.');
                end
            end
        end

        % activate core loss calculation for specfied mesh zones
        function activateCoreLossCalculation(obj, varargin)

            % loop over inputs
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    mzName = obj.m.checkMeshZoneExistence(varargin{i});
                    obj.m.mzs.(mzName).props.isCoreLossActivated = true;
                elseif isvector(varargin{i}) && isstring(varargin{i})
                    for j = 1:numel(varargin{i})
                        mzName = char(varargin{i}(j));
                        mzName = obj.m.checkMeshZoneExistence(mzName);
                        obj.m.mzs.(mzName).props.isCoreLossActivated = true;
                    end
                else
                    error('The input class type must be <char> or <string>.');
                end
            end

        end

        % deactivate core loss calculation for specfied mesh zones
        function deactivateCoreLossCalculation(obj, varargin)

            % loop over inputs
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    mzName = obj.m.checkMeshZoneExistence(varargin{i});
                    obj.m.mzs.(mzName).props.isCoreLossActivated = false;
                elseif isvector(varargin{i}) && isstring(varargin{i})
                    for j = 1:numel(varargin{i})
                        mzName = char(varargin{i}(j));
                        mzName = obj.m.checkMeshZoneExistence(mzName);
                        obj.m.mzs.(mzName).props.isCoreLossActivated = false;
                    end
                else
                    error('The input class type must be <char> or <string>.');
                end
            end

        end

        function y = evalCoreLossModel1(obj, Kh, Ke)
            y = emdlab_flib_calculateIronLossesModel1(obj, Kh, 1, 2, Ke);
        end

        %% Boundary Conditions
        function setAzBC(obj, index, value, varargin)
            obj.bcs.setDirichlet(index, value, varargin{:});
        end

        function clearAllBCs(obj)
            obj.bcs.clearAllBCs;
        end

        function setOddPeriodicBC(obj, varargin)
            obj.bcs.setOddPeriodic(varargin{:});
        end

        function setEvenPeriodicBC(obj, varargin)
            obj.bcs.setEvenPeriodic(varargin{:});
        end
        
        %% Visualaizing Functions
        function varargout = showCoil(obj, coilName)

            coilName = obj.checkCoilExistence(coilName);
            f = emdlab_r2d_mesh;
            ax = axes(f);
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
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = plotMvecOnCenterOfElements(obj)

            f = figure;
            ax = axes(f);
            f.Name = 'Magnetization vectors';
            ax.NextPlot = 'add';
            c = obj.m.getCenterOfElements;
            color = zeros(obj.m.Ne, 3);
            mzNames = fieldnames(obj.m.mzs);

            for i = 1:obj.m.Nmzs
                mzptr = obj.m.mzs.(mzNames{i});
                color(obj.m.ezi(:, mzptr.zi), :) = repmat(mzptr.color, mzptr.Ne, 1);
            end

            patch(ax, 'faces', obj.m.cl(:, 1:3), 'vertices', obj.m.nodes, ...
                'facecolor', 'flat', 'FaceVertexCData', color, 'EdgeColor', 'w', 'facealpha', 1);
            quiver(ax, c(:, 1), c(:, 2), obj.edata.MagnetizationX', obj.edata.MagnetizationY', 'color', 'k');
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', [130,130,130]/255);
            axis(ax, 'off', 'equal');
            zoom on;
            set(f, 'Visible', 'on');
            set(ax, 'clipping', 'off');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = showAllCoilArms(obj)

            f = emdlab_r2d_mesh;
            ax = axes(f);
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
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        %% Solver: Time Settings
        function setSimulationTime(obj, timeStep, totalTime)
            obj.solverSettings.timeStep = timeStep;
            obj.solverSettings.Nsteps = ceil(totalTime / timeStep);
        end

        %% Solver: Auxiliary Settings
        function saveBe(obj, value)
            obj.solverSettings.saveBe = value;
        end

        %% Postproccessing: Evaluations
        %% Postproccessing: plots
        function plotCoilCurrents(obj)

            figure;
            hold on; box on;

            coilNames = fieldnames(obj.coils);

            legNames = {};
            for i = 1:obj.Ncoils
                cptr = obj.coils.(coilNames{i});
                if ~cptr.isCageMember
                    plot(obj.simTime, cptr.current);
                    legNames{end+1} = coilNames{i};
                end
            end

            legend(legNames);
            xlabel('Time')
            ylabel('Coil currents');
            set(gca, 'xlim', [0,obj.simTime(end)]);

        end

        function plotCoilVoltages(obj)

            figure;
            hold on; box on;

            coilNames = fieldnames(obj.coils);

            legNames = {};
            for i = 1:obj.Ncoils
                cptr = obj.coils.(coilNames{i});
                if ~cptr.isCageMember
                    plot(obj.simTime, cptr.voltage);
                    legNames{end+1} = coilNames{i};
                end
            end

            legend(legNames);
            xlabel('Time')
            ylabel('Coil currents');
            set(gca, 'xlim', [0,obj.simTime(end)]);

        end

        function plotCageCurrents(obj, cageName)

            cageName = obj.checkCageExistence(cageName);

            figure;
            hold on; box on;

            coilNames = obj.cages.(cageName).coilArms;

            for i = 1:numel(coilNames)
                cptr = obj.coils.(coilNames(i));
                plot(obj.simTime, cptr.current);
            end

            legend(coilNames);
            xlabel('Time')
            ylabel('Coil currents');
            set(gca, 'xlim', [0,obj.simTime(end)]);

        end

        function plotCoilFluxLinkages(obj)

            figure;
            hold on; box on;

            coilNames = fieldnames(obj.coils);

            legNames = {};
            for i = 1:obj.Ncoils
                cptr = obj.coils.(coilNames{i});
                if ~cptr.isCageMember
                    plot(obj.simTime, cptr.fluxLinkage);
                    legNames{end+1} = coilNames{i};
                end
            end

            legend(legNames);
            xlabel('Time')
            ylabel('Coil flux linkage');
            set(gca, 'xlim', [0,obj.simTime(end)]);
            
        end

        function plotCoilInducedVoltages(obj)

            figure;
            hold on; box on;

            coilNames = fieldnames(obj.coils);

            legNames = {};
            for i = 1:obj.Ncoils
                cptr = obj.coils.(coilNames{i});
                if ~cptr.isCageMember
                    plot(obj.simTime, cptr.inducedVoltage);
                    legNames{end+1} = coilNames{i};
                end
            end

            legend(legNames);
            xlabel('Time')
            ylabel('Coil induced voltage');
            set(gca, 'xlim', [0,obj.simTime(end)]);
            
        end

        function plotWindingsVoltage(obj)
            simTime = 0:obj.solverSettings.timeStep:obj.solverSettings.timeStep * (obj.solverSettings.Nsteps - 1);
            wNames = fieldnames(obj.windings);
            hold all;

            for i = 1:numel(wNames)
                wptr = obj.windings.(wNames{i});
                plot(simTime, obj.results.windingsVoltage(wptr.wi, :));
            end

            legend(wNames);
        end

        function clearAllResults(obj)

            % result names
            rNames = fieldnames(obj.results);

            for i = 1:numel(rNames)
                obj.results.(rNames{i}) = [];
            end

            obj.isResultsValid = false;

        end

        % Evaluation of B [tesla] on gaussian and mesh points of each element
        % 1 gaussian point -> (1/3,1/3)
        % 3 mesh points at corners
        function evalBe(obj)

            [obj.results.Bxg, obj.results.Byg, obj.results.Bxn, obj.results.Byn] = ...
                emdlab_m2d_tl3_evalB(obj.m.cl, obj.results.A, obj.m.JIT);
            obj.results.Bxg = obj.results.Bxg * (obj.units.k_magneticVectorPotential / obj.units.k_length);
            obj.results.Byg = obj.results.Byg * (obj.units.k_magneticVectorPotential / obj.units.k_length);
            obj.results.Bxn = obj.results.Bxn * (obj.units.k_magneticVectorPotential / obj.units.k_length);
            obj.results.Byn = obj.results.Byn * (obj.units.k_magneticVectorPotential / obj.units.k_length);

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

                [obj.results.BxnSmooth(:,eziptr), obj.results.BynSmooth(:,eziptr)] = ...
                    emdlab_m2d_tl3_evalBnSmooth(obj.m.cl(eziptr,:), obj.results.Bxn(:,eziptr), obj.results.Byn(:,eziptr), obj.m.gea(eziptr), obj.m.Nn);

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

                [obj.results.HxnSmooth(:,eziptr), obj.results.HynSmooth(:,eziptr)] = ...
                    emdlab_m2d_tl3_evalBnSmooth(obj.m.cl(eziptr,:), obj.results.Hxn(:,eziptr), obj.results.Hyn(:,eziptr), obj.m.gea(eziptr), obj.m.Nn);

            end

        end

        % calculate total stored energy and co-energy
        function [ye, yc] = evalTotalEnergyCoenergy(obj)

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

        % calculate mesh zone flux linkage
        function y = evalMeshZoneFluxLinkage(obj, mzName)

            mzName = obj.m.checkMeshZoneExistence(mzName);
            y = obj.m.mzs.(mzName).getQ * obj.results.A(obj.m.mzs.(mzName).l2g) * obj.getDepth;

        end

        % calculate winding flux linkage
        function y = evalWindingFluxLinkage(obj, windingName)

            windingName = obj.checkWindingExistence(windingName);
            mptr = obj.exmtcs.(windingName);
            y = 0;

            for i = 1:mptr.Nmzs
                cptr = obj.coils.(mptr.mzsName{i});
                y = y + cptr.sign * cptr.turns * obj.evalMeshZoneFluxLinkage(mptr.mzsName{i});
            end

            y = y / mptr.np;

        end

        function y = evalTorqueByMST(obj, xc, yc, r, N)

            if nargin < 5
                N = 10001;
            end
            [br, bt, t] = obj.getBrBtOnCircle(xc, yc, r, N);
            y = trapz(t, br.*bt) * r^2 * obj.units.k_length^2 * obj.getDepth / (4 * pi * 1e-7);

        end

        function y = evalTorqueByMST3(obj, xc, yc, r, gap, N)

            if nargin < 6, N = 10001; end
            y = (obj.evalTorqueByMST(xc,yc,r-gap/4,N) + obj.evalTorqueByMST(xc,yc,r,N) + obj.evalTorqueByMST(xc,yc,r+gap/4,N))/3;

        end

        function torque = evalTorqueBySurfaceMST(obj, varargin)

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

        function y = evalTorqueByArkkio(obj, mzName, GapLength, xc, yc)

            if nargin < 4
                xc = 0;
                yc = 0;
            end

            obj.evalBe;
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            b = [obj.results.Bxg(obj.m.ezi(:, mzptr.zi)); obj.results.Byg(obj.m.ezi(:, mzptr.zi))]';
            p = mzptr.getCenterOfElements;
            p(:,1) =  p(:,1)-xc;
            p(:,2) =  p(:,2)-yc;
            b = ext_xy2rt(p, b);
            r = sqrt(sum(p.^2, 2));
            y = mzptr.getAreaOfElements' * (b(:, 1) .* b(:, 2) .* r);
            y = y * obj.getDepth * obj.units.k_length^2 / GapLength / (4 * pi * 1e-7);

        end

        function y = evalTorqueByArkkioSmooth(obj, mzName, GapLength, xc, yc)

            if nargin < 4
                xc = 0;
                yc = 0;
            end

            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            bx = sum(obj.results.BxnSmooth(:,obj.m.ezi(:, mzptr.zi)),1)'/3;
            by = sum(obj.results.BynSmooth(:,obj.m.ezi(:, mzptr.zi)),1)'/3;
            b = [bx,by];
            p = mzptr.getCenterOfElements;
            p(:,1) =  p(:,1)-xc;
            p(:,2) =  p(:,2)-yc;
            b = ext_xy2rt(p, b);
            r = sqrt(sum(p.^2, 2));
            y = mzptr.getAreaOfElements' * (b(:, 1) .* b(:, 2) .* r);
            y = y * obj.getDepth * obj.units.k_length^2 / GapLength / (4 * pi * 1e-7);

        end

        function varargout = plotBmag(obj, varargin)

            % amplitude of the B at mesh points
            ti = obj.m.getti(varargin{:});
            ampB = sqrt(obj.results.Bxn(:,ti).^2 + obj.results.Byn(:,ti).^2);

            % plot using patch function
            f = figure;
            ax = axes(f);
            f.Name = 'Magnetic flux density amplitude';

            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(ti, [1,2,3]))', ...
                'CData', ampB, 'FaceColor', 'interp', ...
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

        function varargout = plotCoreLossDensity(obj, mzName, ampB)

            % amplitude of the B at mesh points
            ti = obj.m.getti(mzName);
            mzptr = obj.m.mzs.(mzName);
            eziptr = obj.m.ezi(:,mzptr.zi);
            ampB = repmat(ampB,3,1);
            [ampB, ~] = emdlab_m2d_tl3_evalBnSmooth(obj.m.cl(eziptr,:), ampB, ampB, obj.m.gea(eziptr), obj.m.Nn);

            % plot using patch function
            f = figure;
            ax = axes(f);
            f.Name = 'Core loss density';

            xdata = obj.m.nodes(:,1);
            ydata = obj.m.nodes(:,2);
            patch('XData', xdata(obj.m.cl(ti, [1,2,3]))', 'YData', ydata(obj.m.cl(ti, [1,2,3]))', ...
                'CData', ampB, 'FaceColor', 'interp', ...
                'EdgeColor', 'none', 'parent', ax);

            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k');

            colormap(jet(15));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Loss Density [W/m^3]';

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
            clim([0,1.9]);

            cRange = linspace(min(obj.results.A), max(obj.results.A), N+2);
            c = tmzpc_contour_tl3(obj.m.cl, obj.m.nodes, obj.results.A, cRange(2:end-1));
            t = 1:size(c, 1);
            t = reshape(t, 2, [])';
            patch(ax, 'faces', t, 'vertices', c, 'edgecolor', 'k');

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
        
        function varargout = plotFQvecOnCenterOfElements(obj, FQname)

            f = figure;
            ax = axes(f);            
            ax.NextPlot = 'add';

            switch FQname
                case 'B'
                    obj.evalBn;
                    f.Name = 'Magnetic flux density on center of mesh elements';
                    FQx = obj.results.Bxg';
                    FQy = obj.results.Byg';

                case 'M'
                    f.Name = 'Magnetization vectors on center of mesh elements';
                    FQx = obj.edata.MagnetizationX';
                    FQy = obj.edata.MagnetizationY';
                    
                case 'H'
                    obj.evalHn;
                    f.Name = 'Magnetic field intensity on center of mesh elements';
                    FQx = obj.results.Hxg';
                    FQy = obj.results.Hyg';

                otherwise
                    error('Wrong field quantity.');
            end
            
            % get center of elements
            c = obj.m.getCenterOfElements;
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
            quiver(ax, c(:, 1), c(:, 2), FQx, FQy, 'color', 'k');

            % plot wireframe geometry
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', [130,130,130]/255);

            axis(ax, 'off', 'equal');
            zoom on;
            set(f, 'Visible', 'on');
            set(ax, 'clipping', 'off');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function plotBvecOnCenterOfElements(obj)
            obj.plotFQvecOnCenterOfElements('B');
        end

        function plotHvecOnCenterOfElements(obj)
            obj.plotFQvecOnCenterOfElements('H');
        end

%         function plotMvecOnCenterOfElements(obj)
%             obj.plotFQvecOnCenterOfElements('M');
%         end

        %% post-proccessing: interpolation functions
        % get triangle index of points
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
            ti = pointLocation(triangulation(obj.m.cl, obj.m.nodes), [x,y]);

        end

        % get the value of Az on points
        function Az = getAzOnPoints(obj, x, y)

            % find index of triangle containing the (x,y) point
            ti = obj.getPointsLocation(x, y);

            % care points out of the mesh
            index = isnan(ti);
            ti(index) = 1;

            % interpolate the value of A on points
            Az = emdlab_m2d_tl3_interpA(obj.m.cl, obj.m.nodes, obj.results.A(obj.m.cl'), obj.m.JIT, ti, x, y);
           
            Az(index) = NaN;

            if iscolumn(Az), Az = Az'; end

        end

        % get x- and y-components of a field quantity: B or H
        function [FQx, FQy] = getFQxFQyOnPoints(obj, x, y, fieldQuantity)

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

        % get the value of Bx and By on points that are on a segment
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

        % get the value of Bx and By on points that are on a segment
        function [Hx, Hy, x, y] = getHxHyOnSegment(obj, x1, y1, x2, y2, N)

            % set default value of N
            if nargin < 6, N = 1000; end
            [Hx, Hy, x, y] = obj.getFQxFQyOnSegment(x1, y1, x2, y2, N, 'HSmooth');

        end

        % get the value of Bx and By on points that are on a circle
        function [FQx, FQy, x, y] = getFQxFQyOnCircle(obj, xc, yc, r, N, fieldQuantity)

            % sample points
            t = linspace(0,2*pi,N);
            x = xc + r * cos(t);
            y = yc + r * sin(t);
            [FQx, FQy] = obj.getFQxFQyOnPoints(x, y, fieldQuantity);

        end

        % get the value of Bx and By on points that are on a circle
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

        % get the value of Bx and By on points that are on a circle
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

        % get the value of Br and Bt on points that are on a circle
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

        % get the value of Br and Bt on points that are on a circle
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

        function [Bx,By] = getMeshZoneBxgByg(obj, mzName)
            Bx = obj.results.Bxg(obj.m.ezi(:,obj.m.mzs.(mzName).zi));
            By = obj.results.Byg(obj.m.ezi(:,obj.m.mzs.(mzName).zi));
        end

        function Az = getMeshZoneAz(obj, mzName)
            Az = obj.results.A(obj.m.mzs.(mzName).l2g);
        end

    end

    %% Getters
    methods

        function y = getDepth(obj)
            y = obj.depth * obj.units.getQuantityScaler('length');
        end

    end

    %% Internal Methods
    methods (Access = protected)

        function makeFalse_isElementDataAssigned(obj)
            obj.isElementDataAssigned = false;
        end

    end

end

% EMDLAB: Electrical Machines Design Laboratory
% A 3-dimensional thermal-static solver based on thermal network theory

classdef emdlab_solvers_ts3d_tn < handle

    properties (SetAccess = protected)

        % solver mesh
        m (1,1);

        % internal & boundary conditions
        excitations (1,1) struct;

        % contacts
        contacts (1,1) struct;

        % elements data
        edata (1,1) struct;

        % results
        results (1,1) struct;

        % matrix containing source terms due to aisotropy
        kxyzSV (:,:) double;

        % initial temperature assumption for iterative solver
        initialTemperature (1,1) double = 25;

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

        meshType (1,:) char = '';
        NR (1,1) double; % number of resistances per element
        Nfn (1,1) double; % maximum number of nodes of one mesh facet

    end

    properties (Dependent = true)

        % number of excitations
        Nex (1,1) double;

    end

    methods
        %% Constructor and Destructor
        function obj = emdlab_solvers_ts3d_tn(m)

            % mesh pointer
            switch class(m)
                case 'emdlab_m3d_thmdb'
                    obj.meshType = 'thm';
                    obj.NR = 4;
                    obj.Nfn = 3;
                case 'emdlab_m3d_hhmdb'
                    obj.meshType = 'hhm';
                    obj.NR = 6;
                    obj.Nfn = 4;
                case 'emdlab_m3d_pmdb'
                    obj.meshType = 'pm';
                    obj.NR = 5;
                    obj.Nfn = 4;
                otherwise
                    error('Mesh type is not supported.');
            end
            m.ggmesh; % generate global mesh
            obj.m = m;
            if isa(m, 'emdlab_m3d_thmdb'); m.evalJIT; end

            % set default properties of mesh zones
            for mzName = obj.m.getMeshZoneNames
                obj.setdp(mzName);
            end

            % default values
            obj.units = emdlab_phy_units;

        end

        function delete(obj)
            delete(obj.m);
            delete(obj.units);
        end

        %% Unit Manager
        function setLengthUnit(obj, unitValue)
            obj.units.setQuantityUnit('length', unitValue);
        end

        function setUnit(obj, varargin)
            obj.units.setQuantityUnit(varargin{:});
        end

        %% Design Properties
        function setInitialTemperatureOfAllMeshZones(obj, value)
            obj.initialTemperature = value;
        end

        %% solver properties for mesh zones

        function setdp(obj, mzName)
            % set default properties of a mesh zone

            obj.m.mzs.(mzName).props.initialTemperature = 25;

        end

        function rotateMeshZone(obj, mzName, varargin)

            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.m.rotateMeshZone(mzName, varargin{:});

            if obj.m.mzs.(mzName).props.isMagnetized
                obj.m.mzs.(mzName).props.magnetization.rotate(varargin{:});
            end

        end

        %% Solver Functions
        function setSolverIndex(obj, value)
            % solver index used for parallel computing
            obj.solverIndex = value;
        end

        function assignElementData(obj)

            % check states
            if obj.isElementDataAssigned, return; end

            % allocation of memory
            obj.edata.ThermalConductivity = zeros(3, obj.m.Ne);

            % getting mesh zones
            mzNames = obj.m.getMeshZoneNames;

            % check if all materials are linear
            obj.edata.areAllTemperatureIndependent = true;
            obj.edata.areAllIsotropic = true;
            obj.edata.areAllHomogeneous = true;

            % loop over mesh zones
            for i = 1:obj.m.Nmzs

                % get pointer to mesh zone
                mzptr = obj.m.mzs.(mzNames{i});

                % assigning thermal conductivities
                if obj.m.mts.(mzptr.material).ThermalConductivity.isTemperatureDependent
                    if obj.m.mts.(mzptr.material).ThermalConductivity.isIsotropic
                        obj.edata.ThermalConductivity(:,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ThermalConductivity.value;
                    else
                        obj.edata.ThermalConductivity(:,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ThermalConductivity.value(:);
                    end
                else
                    if obj.m.mts.(mzptr.material).ThermalConductivity.isIsotropic
                        obj.edata.ThermalConductivity(:,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ThermalConductivity.value;
                    else
                        obj.edata.ThermalConductivity(1,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ThermalConductivity.value(1);
                        obj.edata.ThermalConductivity(2,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ThermalConductivity.value(2);
                        obj.edata.ThermalConductivity(3,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ThermalConductivity.value(3);
                        obj.edata.areAllIsotropic = false;
                    end
                end

            end

            % change states
            obj.isElementDataAssigned = true;

        end

        function solve(obj)

            % assign element data
            obj.assignElementData;

            % calculate resistances
            resistances = zeros(obj.m.Ne,obj.NR);
            um = obj.units.k_length;

            % loop over mesh zones
            mzNames = obj.m.getMeshZoneNames;
            ux = zeros(obj.m.Ne,3);
            uy = zeros(obj.m.Ne,3);
            uz = zeros(obj.m.Ne,3);
            for mz = mzNames
                mzptr = obj.m.mzs.(mz);
                ux(obj.m.ezi(:,mzptr.zi),:) = repmat(obj.m.cs.(mzptr.orientation).ux,mzptr.Ne,1);
                uy(obj.m.ezi(:,mzptr.zi),:) = repmat(obj.m.cs.(mzptr.orientation).uy,mzptr.Ne,1);
                uz(obj.m.ezi(:,mzptr.zi),:) = repmat(obj.m.cs.(mzptr.orientation).uz,mzptr.Ne,1);
            end

            % store source term due to anisotropy
            obj.kxyzSV = zeros(obj.m.Ne,obj.NR);

            % distance modifier
            elm = zeros(obj.m.Ne,obj.NR);

            % loop over elements to calculate resistances of each element
            for i = 1:obj.m.Ne

                kx = obj.edata.ThermalConductivity(1,i);
                ky = obj.edata.ThermalConductivity(2,i);
                kz = obj.edata.ThermalConductivity(3,i);

                for j = 1:obj.m.elements(i,end)

                    % direction vector from cell center i to cell center j
                    if obj.m.nbs(i,j)
                        uij = obj.m.elementCenter(obj.m.nbs(i,j), :) - obj.m.elementCenter(i,:);
                        die = 0.5 * norm(uij);
                    else
                        uij = obj.m.facetCenter(abs(obj.m.elements(i,j)),:) - obj.m.elementCenter(i,:);
                        die = norm(uij);
                    end
                    uij = uij / norm(uij);

                    kij = kx*(ux(i,:)*uij').^2 + ky*(uy(i,:)*uij').^2 + kz*(uz(i,:)*uij').^2;

                    elm(i,j) = abs(uij * obj.m.facetNormal(abs(obj.m.elements(i,j)),:)');

                    % calculation of conductances
                    resistances(i,j) = die/(kij * elm(i,j) * obj.m.fa(i,j) * um);

                end

            end

            % initialize conductance matrix and source vector
            [indexI,indexJ,value,sourceVector] = buildGMatrixSVector(obj, resistances);

            % apply excitation conditions
            exNames = obj.getExcitationNames;
            for i = 1:numel(exNames)

                exptr = obj.excitations.(exNames(i));
                switch exptr.type
                    case 'fixed-temperature'
                        [value, sourceVector] = applyFixedTemperatureBC(obj, exptr.idx, value, sourceVector, resistances, exptr.value);
                    case 'convection'
                        [value, sourceVector] = applyConvectionBC(obj, exptr.idx, value, sourceVector, resistances, exptr.Tinf, exptr.hValue);
                    case 'radiation'
                        [value, sourceVector] = applyConvectionBC(obj, exptr.idx, value, sourceVector, resistances, exptr.Tinf, exptr.hValue);
                    case 'heat-flux'
                        sourceVector = applyHeatFluxBC(obj, exptr.idx, sourceVector, exptr.value);
                    case 'internal-heat-source'
                        sourceVector = applyIHG(obj, exptr.mzName, sourceVector, exptr.value);
                end

            end

            % apply contacts
            cNames = obj.getContactsNames;
            for i = 1:numel(cNames)
                value = obj.applyContact(obj.contacts.(cNames(i)).idx, value, obj.contacts.(cNames(i)).hValue);
            end

            % construct global matrices and solve
            G_maxtrix = sparse(indexI, indexJ, value);
            obj.results.T = G_maxtrix\sourceVector;
            obj.evalTn;

            obj.solverHistory.relativeError = [];

            %             if ~obj.edata.areAllTemperatureIndependent || ~obj.edata.areAllIsotropic || ~obj.edata.areAllHomogeneous
            %
            %                 releativeError = 1e-3;
            %                 maxIteration = 200;
            %                 iter = 0;
            %                 err = inf;
            %
            %                 % iterative loop for solver
            %                 while iter<maxIteration && err>releativeError
            %
            %                     sourceVectorU = sourceVector;
            %
            %                     % handle anisotropy
            %                     if ~obj.edata.areAllIsotropic
            %                         for i = 1:obj.m.Ne
            %                             Tij = obj.results.Tn(obj.m.cl(i,[1:obj.NR,1]));
            %                             for j = 1:obj.NR
            %                                 if ~obj.m.bedges(abs(obj.m.elements(i,j)))
            %                                     sourceVectorU(i) = sourceVectorU(i) + z * ...
            %                                         obj.kxySV(i,j) * (Tij(j) - Tij(j+1)) * elm(i,j);
            %                                 end
            %                             end
            %                         end
            %                     end
            %
            %                     Told = obj.results.T;
            %
            %                     % solver iteration
            %                     obj.results.T = G_maxtrix\sourceVectorU;
            %                     obj.evalTn;
            %
            %                     iter = iter + 1;
            %                     err = norm(obj.results.T-Told,2)/norm(obj.results.T,2);
            %                     obj.solverHistory.relativeError(end+1) = err;
            %
            %                     fprintf('Iteration #%03d, Relative Error = %.2e\n', iter, err);
            %
            %                 end
            %
            %             end

            % calculate cross heat at facets
            Nfn_ = obj.Nfn;
            Nfacets = size(obj.m.facets,1);
            obj.results.qf = zeros(Nfacets,1);
            %             for k = 1:Nfacets
            %                 if ~obj.m.bfacets(k)
            %                     % index of left element
            %                     i = obj.m.facets(k,Nfn_+3);
            %                     % index of right element
            %                     j = obj.m.facets(k,Nfn_+5);
            %                     % net heat passing k-th internal edge
            %                     obj.results.qf(k) = (obj.results.T(j) - obj.results.T(i)) * G_maxtrix(i,j);
            %                 end
            %             end

            emdlab_mex_solvers_ts3d_calcqfif(obj.m.facets, obj.m.bfacets, obj.results.T, G_maxtrix, obj.results.qf, Nfn_);

            for k = 1:numel(exNames)

                exptr = obj.excitations.(exNames(k));
                switch exptr.type

                    case 'fixed-temperature'
                        if isa(exptr.value, 'function_handle')
                            Tb = exptr.value(obj.m.facetCenter(exptr.idx,1),obj.m.facetCenter(exptr.idx,2),obj.m.facetCenter(exptr.idx,3));
                        else
                            Tb = exptr.value*ones(length(exptr.idx),1);
                        end

                        emdlab_mex_solvers_ts3d_calcqfft(resistances, obj.m.facets, exptr.idx, obj.results.T, Tb, obj.results.qf, Nfn_);

                        %                         j = 0;
                        %                         for i = exptr.idx(:)'
                        %                             j = j+1;
                        %
                        %                             if obj.m.facets(i,Nfn_+2)
                        %                                 obj.results.qf(i) = (Tb(j) - obj.results.T(obj.m.facets(i,Nfn_+5))) / ...
                        %                                     resistances(obj.m.facets(i,Nfn_+5), obj.m.facets(i,Nfn_+6));
                        %                             else
                        %                                 obj.results.qf(i) = (obj.results.T(obj.m.facets(i,Nfn_+3)) - Tb(j)) / ...
                        %                                     resistances(obj.m.facets(i,Nfn_+3), obj.m.facets(i,Nfn_+4));
                        %                             end
                        %                         end

                    case 'convection'
                        if isa(exptr.hValue, 'function_handle')
                            hValue = exptr.hValue(edgeCenter(exptr.idx,1),edgeCenter(exptr.idx,2));
                        else
                            hValue = exptr.hValue*ones(length(exptr.idx),1);
                        end

                        j = 0;
                        for i = exptr.idx(:)'
                            j = j+1;

                            if obj.m.facets(i,Nfn_+2)
                                obj.results.qf(i) = (exptr.Tinf - obj.results.T(obj.m.facets(i,Nfn_+5))) / ...
                                    (resistances(obj.m.facets(i,Nfn_+5), obj.m.facets(i,Nfn_+6)) + 1/(hValue(j)*obj.m.facetArea(i)*obj.units.k_length^2));
                            else
                                obj.results.qf(i) = (obj.results.T(obj.m.facets(i,Nfn_+3)) - exptr.Tinf) / ...
                                    (resistances(obj.m.facets(i,Nfn_+3), obj.m.facets(i,Nfn_+4)) + 1/(hValue(j)*obj.m.facetArea(i)*obj.units.k_length^2));
                            end
                        end

                    case 'radiation'
                        [value, sourceVector] = applyConvectionBC(obj, exptr.idx, value, sourceVector, resistances, exptr.Tinf, exptr.hValue);

                    case 'heat-flux'
                        if isa(exptr.value, 'function_handle')
                            qb = exptr.value(obj.m.facetCenter(exptr.idx,1),obj.m.facetCenter(exptr.idx,2),obj.m.facetCenter(exptr.idx,3));
                        else
                            qb = exptr.value*ones(length(exptr.idx),1);
                        end

                        emdlab_mex_solvers_ts3d_calcqfhf(obj.m.facets, exptr.idx, qb, obj.results.qf, obj.m.facetArea, Nfn_, obj.units.k_length^2);

                        %                         j = 0;
                        %                         for i = exptr.idx(:)'
                        %                             j = j+1;
                        %
                        %                             if obj.m.facets(i,Nfn_+2)
                        %                                 obj.results.qf(i) = qb(j)*obj.m.facetArea(i)*obj.units.k_length^2;
                        %                             else
                        %                                 obj.results.qf(i) = -qb(j)*obj.m.facetArea(i)*obj.units.k_length^2;
                        %                             end
                        %                         end

                end

            end

            obj.evalTSmooth;
        end

        function y = calculateNetHeatCrossingBoundaryFacets(obj, idx)
            y = 0;
            Nfn_ = obj.Nfn;
            for i = idx(:)'
                if obj.m.facets(i,Nfn_+2)
                    y = y + obj.results.qf(i);
                else
                    y = y - obj.results.qf(i);
                end
            end
        end

        function evalTn(obj)
            obj.results.Tn = emdlab_mex_solvers_ts3d_evalTn(obj.results.T, obj.m.elementCenter, obj.m.nodes, obj.m.cl);
        end

        function evalTn_mfile(obj)

            obj.results.Tn = zeros(obj.m.Nn,1);
            ec = obj.m.elementCenter;
            denum = zeros(obj.m.Nn,1);

            % loop over elements
            for i = 1:obj.m.Ne
                for j = 1:obj.m.cl(i,end)
                    pIndex = obj.m.cl(i,j);
                    di = norm(ec(i,:) - obj.m.nodes(pIndex,:));
                    obj.results.Tn(pIndex) = obj.results.Tn(pIndex) + obj.results.T(i)/di;
                    denum(pIndex) = denum(pIndex) + 1/di;
                end
            end

            obj.results.Tn = obj.results.Tn ./ denum;

            %             % apply excitation conditions
            %             exNames = obj.getExcitationNames;
            %             for i = 1:numel(exNames)
            %                 exptr = obj.excitations.(exNames(i));
            %                 if strcmpi(exptr.type,'fixed-temperature')
            %                     for j = reshape(exptr.idx, 1, [])
            %                         if isa(exptr.value, 'function_handle')
            %                             pts = obj.m.nodes(obj.m.edges(j,1:2),:);
            %                             obj.results.Tn(obj.m.edges(j,1)) = exptr.value(pts(1,1),pts(1,2));
            %                             obj.results.Tn(obj.m.edges(j,2)) = exptr.value(pts(2,1),pts(2,2));
            %                         else
            %                             obj.results.Tn(obj.m.edges(j,1:2)) = exptr.value;
            %                         end
            %                     end
            %                 end
            %             end

        end

        function evalTSmooth(obj)

            obj.results.Tsmooth = obj.results.Tn;

            % apply excitation conditions
            exNames = obj.getExcitationNames;
            for i = 1:numel(exNames)
                exptr = obj.excitations.(exNames(i));
                if strcmpi(exptr.type,'fixed-temperature')
                    for j = exptr.idx(:)'
                        if isa(exptr.value, 'function_handle')
                            pts = obj.m.nodes(obj.m.facets(j,1:obj.m.facets(j,end)),:);
                            for k = 1:obj.m.facets(j,end)
                                obj.results.Tsmooth(obj.m.facets(j,k)) = exptr.value(pts(k,1),pts(k,2),pts(k,3));
                            end
                        else
                            obj.results.Tsmooth(obj.m.facets(j,1:obj.m.facets(j,end))) = exptr.value;
                        end
                    end
                end
            end

        end

        %% Excitation Definitions
        function addFixedTemperatureBC(obj, exName, idx, value)

            exName = obj.checkExcitationNonExistence(exName);
            obj.excitations.(exName).type = 'fixed-temperature';
            obj.excitations.(exName).idx = idx;
            obj.excitations.(exName).value = value;

        end

        function addConvectionBC(obj, exName, idx, hValue, Tinf)

            exName = obj.checkExcitationNonExistence(exName);
            obj.excitations.(exName).type = 'convection';
            obj.excitations.(exName).idx = idx;
            obj.excitations.(exName).hValue = hValue;
            obj.excitations.(exName).Tinf = Tinf;

        end

        function addRadiationBC(obj, exName, idx, eValue, Tsurr)

            exName = obj.checkExcitationNonExistence(exName);
            obj.excitations.(exName).type = 'radiation';
            obj.excitations.(exName).idx = idx;
            obj.excitations.(exName).eValue = eValue;
            obj.excitations.(exName).Tsurr = Tsurr;

        end

        function addHeatFluxBC(obj, exName, idx, value)

            exName = obj.checkExcitationNonExistence(exName);
            obj.excitations.(exName).type = 'heat-flux';
            obj.excitations.(exName).idx = idx;
            obj.excitations.(exName).value = value;

        end

        function addInternalHeatSource(obj, exName, mzName, value, unit)

            exName = obj.checkExcitationNonExistence(exName);
            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.excitations.(exName).type = 'internal-heat-source';
            obj.excitations.(exName).mzName = mzName;


%             z = obj.getDepth;
            um = obj.units.k_length;
%             Area = sum(obj.m.gea(obj.m.ezi(:,obj.m.mzs.(mzName).zi))) * um^2;

            % set default unit
            if nargin < 5, unit = 'W'; end

            unit = lower(erase(unit, ' '));
            switch unit

                case 'w'
                    value = value/(obj.m.mzs.(mzName).volume * um^3);

                case {'w/m^3', 'w/m3', 'w/(m^3)'}
                otherwise
                    error('Unsupported unit type.');

            end

            obj.excitations.(exName).value = value;

        end

        function modifyFixedTemperatureBC(obj, exName, varargin)
        % 1. Verify that the boundary condition actually exists first
        if ~isfield(obj.excitations, exName)
            error('EMDLAB:ThermalSolver:ExcitationNotFound', ...
                'Excitation "%s" does not exist. Use addFixedTemperatureBC first.', exName);
        end
        
        % 2. Ensure we are modifying a fixed-temperature BC
        if ~strcmp(obj.excitations.(exName).type, 'fixed-temperature')
            error('EMDLAB:ThermalSolver:InvalidType', ...
                'Excitation "%s" is not a fixed-temperature boundary condition.', exName);
        end

        % 3. Parse optional inputs to allow updating either 'value', 'idx', or both
        p = inputParser;
        p.addParameter('value', [], @isnumeric);
        p.addParameter('idx', [], @isnumeric);
        p.parse(varargin{:});
        
        % Update only the parameters that were provided
        if ~isempty(p.Results.value)
            obj.excitations.(exName).value = p.Results.value;
        end
        if ~isempty(p.Results.idx)
            obj.excitations.(exName).idx = p.Results.idx;
        end
        end

        function modifyRadiationBC(obj, exName, varargin)
        % 1. Verify that the boundary condition actually exists first
        if ~isfield(obj.excitations, exName)
            error('EMDLAB:ThermalSolver:ExcitationNotFound', ...
                'Excitation "%s" does not exist. Use addRadiationBC first.', exName);
        end
        
        % 2. Ensure we are modifying a radiation BC
        if ~strcmp(obj.excitations.(exName).type, 'radiation')
            error('EMDLAB:ThermalSolver:InvalidType', ...
                'Excitation "%s" is not a radiation boundary condition.', exName);
        end

        % 3. Parse optional inputs to allow updating eValue, Tsurr, idx, or combinations of them
        p = inputParser;
        p.addParameter('eValue', [], @isnumeric);
        p.addParameter('Tsurr', [], @isnumeric);
        p.addParameter('idx', [], @isnumeric);
        p.parse(varargin{:});
        
        % Update only the parameters that were provided
        if ~isempty(p.Results.eValue)
            obj.excitations.(exName).eValue = p.Results.eValue;
        end
        if ~isempty(p.Results.Tsurr)
            obj.excitations.(exName).Tsurr = p.Results.Tsurr;
        end
        if ~isempty(p.Results.idx)
            obj.excitations.(exName).idx = p.Results.idx;
        end
        end

        function modifyConvectionBC(obj, exName, hValue, Tinf)
        % 1. Verify that the boundary condition actually exists first
        if ~isfield(obj.excitations, exName)
            error('EMDLAB:ThermalSolver:ExcitationNotFound', ...
                'Excitation "%s" does not exist. Use addConvectionBC first.', exName);
        end
        

        obj.excitations.(exName).hValue =hValue;

            obj.excitations.(exName).Tinf = Tinf;
        end

        function modifyHeatFluxBC(obj, exName, varargin)
        % 1. Verify that the boundary condition actually exists first
        if ~isfield(obj.excitations, exName)
            error('EMDLAB:ThermalSolver:ExcitationNotFound', ...
                'Excitation "%s" does not exist. Use addHeatFluxBC first.', exName);
        end
        
        % 2. Ensure we are modifying a heat-flux BC
        if ~strcmp(obj.excitations.(exName).type, 'heat-flux')
            error('EMDLAB:ThermalSolver:InvalidType', ...
                'Excitation "%s" is not a heat-flux boundary condition.', exName);
        end

        % 3. Parse optional inputs to allow updating value, idx, or both
        p = inputParser;
        p.addParameter('value', [], @isnumeric);
        p.addParameter('idx', [], @isnumeric);
        p.parse(varargin{:});
        
        % Update only the parameters that were provided
        if ~isempty(p.Results.value)
            obj.excitations.(exName).value = p.Results.value;
        end
        if ~isempty(p.Results.idx)
            obj.excitations.(exName).idx = p.Results.idx;
        end
    end

    function modifyInternalHeatSource(obj, exName, value)
        
        mzName = obj.excitations.(exName).mzName;

        um = obj.units.k_length;
%             Area = sum(obj.m.gea(obj.m.ezi(:,obj.m.mzs.(mzName).zi))) * um^2;

            % set default unit
            if nargin < 4, unit = 'W'; end

            unit = lower(erase(unit, ' '));
            switch unit

                case 'w'
                    value = value/(obj.m.mzs.(mzName).volume * um^3);

                case {'w/m^3', 'w/m3', 'w/(m^3)'}
                otherwise
                    error('Unsupported unit type.');

            end

            obj.excitations.(exName).value = value;

    end

        %% Contacts
        function addContact(obj, cName, mz1Name, mz2Name, hValue)

            cName = obj.checkContactNonExistence(cName);
            mz1Name = obj.m.checkMeshZoneExistence(mz1Name);
            mz2Name = obj.m.checkMeshZoneExistence(mz2Name);
            obj.contacts.(cName).idx = obj.m.getFacetIndicesOnContact(mz1Name,mz2Name);
            obj.contacts.(cName).mz1Name = mz1Name;
            obj.contacts.(cName).mz2Name = mz2Name;
            obj.contacts.(cName).hValue = hValue;

        end

        function area = getContactArea(obj, mz1Name, mz2Name)
            idx = obj.m.getFacetIndicesOnContact(mz1Name, mz2Name);
            area = obj.getDepth * sum(obj.m.edgeLength(idx)) * obj.units.k_length;
        end

        %% Visualization Functions
        function varargout = plotAverageTemperature(obj, N, varargin)

            [f,ax] = emdlab_r3d_geometry(1,0);
            if nargin<2, N=10; end

            switch obj.meshType
                case 'thm'
                    idx = obj.m.facets(obj.m.bfacets,6) + obj.m.facets(obj.m.bfacets,8);
                    patch('faces', obj.m.facets(obj.m.bfacets,1:3), 'Vertices', obj.m.nodes, 'FaceVertexCData',obj.results.T(idx), ...
                        'FaceColor','flat', 'edgecolor', 'w');
                case 'hhm'
                    idx = obj.m.facets(obj.m.bfacets,7) + obj.m.facets(obj.m.bfacets,9);
                    patch('faces', obj.m.facets(obj.m.bfacets,1:4), 'Vertices', obj.m.nodes, 'FaceVertexCData',obj.results.T(idx), ...
                        'FaceColor','flat', 'edgecolor', 'w');
            end

            colormap(jet(N));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Average Temperature [C]';
            climits = clim;
            cb.Ticks = fix(linspace(climits(1), climits(2), 10)*100)/100;
            cb.Ticks(1) = cb.Ticks(1) + 0.01;
            cb.Position = [0.8,0.3,0.05,0.4];

            %             index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            %             patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
            %                 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            %             zoom on;
            %             axis(ax, 'off');
            %             axis(ax, 'equal');
            %             set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotHeatFlux(obj, N, varargin)

            [f,ax] = emdlab_r3d_geometry(1,0);
            if nargin<2, N=10; end

            clr = abs(obj.results.qf./obj.m.facetArea)*(1/obj.units.k_length^2);

            switch obj.meshType
                case 'thm'
                    idx = obj.m.facets(obj.m.bfacets,6) + obj.m.facets(obj.m.bfacets,8);
                    patch('faces', obj.m.facets(obj.m.bfacets,1:3), 'Vertices', obj.m.nodes, 'FaceVertexCData',obj.results.T(idx), ...
                        'FaceColor','flat', 'edgecolor', 'w');
                case {'pm','hhm'}
                    idx = (obj.results.qf ~= 0) & obj.m.bfacets;
                    tmp = obj.m.facets(idx,1:4);
                    tmp(tmp == 0) = nan;
                    patch('faces', tmp, 'Vertices', obj.m.nodes, 'FaceVertexCData',clr(idx), ...
                        'FaceColor','flat', 'edgecolor', 'w');
            end

            colormap(jet(N));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Heat Flux [W/m^2]';
            climits = clim;
            cb.Ticks = linspace(climits(1), climits(2), 10)*100/100;
            ticks = cb.Ticks;
            cb.TickLabels = compose('%0.2e',ticks);
            cb.Position = [0.75,0.3,0.05,0.4];
            
            % set outputs
            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotTemperature(obj, N, varargin)

            % set default number of contours
            if nargin < 2, N=10; end

            % get figure & axis
            [f,ax] = emdlab_r3d_geometryNEW(0,0);
            
            % index of boundary facets for mesh zone
            idx = ((obj.m.facets(:,obj.m.NFN+1) == 0) & (obj.m.facets(:,obj.m.NFN+2) ~= 0)) | ...
                ((obj.m.facets(:,obj.m.NFN+1) ~= 0) & (obj.m.facets(:,obj.m.NFN+2) == 0));

            % facet connectivity list
            fcl = obj.m.facets(idx,1:obj.m.NFN);
            fcl(fcl == 0) = nan;

            % plot mesh zone temperature
            patch('faces', fcl, 'Vertices', obj.m.nodes, 'FaceVertexCData',obj.results.Tsmooth, ...
                        'FaceColor','interp', 'edgecolor', 'none');

            % set color bar
            colormap(jet(N));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Temperature [C]';
            climits = clim;
            cb.Ticks = linspace(climits(1), climits(2), 10)*100/100;
            ticks = cb.Ticks;
            cb.TickLabels = compose('%0.2e',ticks);
            cb.Position = [0.75,0.3,0.05,0.4];

            % set outputs
            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotMeshZoneTemperature(obj, mzName, N)

            % set default number of contours
            if nargin < 3, N=10; end

            % check mesh zone name
            mzName = obj.m.checkMeshZoneExistence(mzName);
            zi = obj.m.mzs.(mzName).zi; % zone index

            % get figure & axis
            [f,ax] = emdlab_r3d_geometryNEW(0,0);
            
            % index of boundary facets for mesh zone
            idx = ((obj.m.facets(:,obj.m.NFN+1) == zi) & (obj.m.facets(:,obj.m.NFN+2) ~= zi)) | ...
                ((obj.m.facets(:,obj.m.NFN+1) ~= zi) & (obj.m.facets(:,obj.m.NFN+2) == zi));

            % facet connectivity list
            fcl = obj.m.facets(idx,1:obj.m.NFN);
            fcl(fcl == 0) = nan;

            % plot mesh zone temperature
            patch('faces', fcl, 'Vertices', obj.m.nodes, 'FaceVertexCData',obj.results.Tsmooth, ...
                        'FaceColor','interp', 'edgecolor', 'none');

            % Temperature of mesh zone nodes -> this is used to set color bar
            Tmz = obj.results.Tsmooth(obj.m.mzs.(mzName).l2g);

            % set color bar
            colormap(jet(N));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Temperature [C]';
            clim([min(Tmz), max(Tmz)]);
            climits = clim;
            cb.Ticks = linspace(climits(1), climits(2), 10)*100/100;
            ticks = cb.Ticks;
            cb.TickLabels = compose('%0.2e',ticks);
            cb.Position = [0.75,0.3,0.05,0.4];

            % set outputs
            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end


        function varargout = plotTemperatureTn(obj, N, varargin)

            [f,ax] = emdlab_r3d_geometry(1,0);
            if nargin<2, N=10; end

            switch obj.meshType
                case 'thm'

                    patch('faces', obj.m.facets(obj.m.bfacets,1:3), 'Vertices', obj.m.nodes, 'FaceVertexCData',obj.results.Tn, ...
                        'FaceColor','interp', 'edgecolor', 'w');
                case {'pm','hhm'}
                    tmp = obj.m.facets(obj.m.bfacets,1:4);
                    tmp(tmp ==0 ) = nan;
                    patch('faces', tmp, 'Vertices', obj.m.nodes, 'FaceVertexCData',obj.results.Tn, ...
                        'FaceColor','interp', 'edgecolor', 'w');
            end

            colormap(jet(N));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Average Temperature [C]';
            climits = clim;
            cb.Ticks = fix(linspace(climits(1), climits(2), 10)*100)/100;
            cb.Ticks(1) = cb.Ticks(1) + 0.01;
            cb.Position = [0.8,0.3,0.05,0.4];

            %             index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            %             patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
            %                 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            %             zoom on;
            %             axis(ax, 'off');
            %             axis(ax, 'equal');
            %             set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function y = getMeshZoneTemperature(obj, mzName)
            idx = obj.m.ezi(:,obj.m.mzs.(mzName).zi);
            y = [min(obj.results.T(idx)), max(obj.results.T(idx)), obj.m.gev(idx)*obj.results.T(idx)/sum(obj.m.gev(idx))];
        end

%         function plotHeatFluxVectors(obj)
%             [f,ax] = obj.m.showwf;
%             q = -obj.results.qf .* obj.m.facetNormal;
%             quiver3(obj.m.facetCenter(:,1),obj.m.facetCenter(:,2),obj.m.facetCenter(:,3), ...
%                 q(:,1), q(:,2), q(:,3), 5, 'Parent', ax);
%             axis off equal
%         end

        function y = getAverageTemperature(obj)
            y = obj.m.gev * obj.results.T / sum(obj.m.gev);
        end

        function y = getMinimumTemperature(obj)
            y = min([obj.results.T; obj.results.Tn; obj.results.Tsmooth]);
        end

        function y = getMaximumTemperature(obj)
            y = max([obj.results.T; obj.results.Tn; obj.results.Tsmooth]);
        end

    end

    methods (Access=private)

        function exName = checkExcitationExistence(obj, exName)

            if ~isfield(obj.excitations, exName)
                throw(MException('', ['Excitation with name [', exName, '] does not exist.']));
            end

        end

        function exName = checkExcitationNonExistence(obj, exName)

            if isfield(obj.excitations, exName)
                throw(MException('', ['Another excitation with name [', exName, '] already exist.']));
            end

        end

        function exNames = getExcitationNames(obj)
            exNames = string(fieldnames(obj.excitations))';
        end

        function cName = checkContactExistence(obj, cName)

            if ~isfield(obj.contacts, cName)
                throw(MException('', ['Contatct with name [', cName, '] does not exist.']));
            end

        end

        function cName = checkContactNonExistence(obj, cName)

            if isfield(obj.contacts, cName)
                throw(MException('', ['Another contacts with name [', cName, '] already exist.']));
            end

        end

        function cNames = getContactsNames(obj)
            cNames = string(fieldnames(obj.contacts))';
        end

        function [indexI,indexJ,value,sourceVector] = buildGMatrixSVector(obj, resistances)

            indexI = zeros(obj.m.Ne,obj.NR+1);
            indexJ = zeros(obj.m.Ne,obj.NR+1);
            value = zeros(obj.m.Ne,obj.NR+1);
            sourceVector = zeros(obj.m.Ne,1);
            Rp = zeros(1,obj.NR);
            idx = 1:obj.NR;

            % loop over elements to construct G_matrix
            for i = 1:obj.m.Ne

                % constructing i-th row
                indexI(i,:) = i;

                % conductance between cell P and j-th nb
                for j = 1:obj.m.elements(i,end)

                    % index of j-th neighbor
                    nbIndex = obj.m.nbs(i,j);

                    if nbIndex
                        Rp(j) = resistances(i,j) + resistances(nbIndex,idx(abs(obj.m.elements(nbIndex,1:obj.NR)) == (abs(obj.m.elements(i,j)))));
                        indexJ(i,j+1) = nbIndex;
                        value(i,j+1) = -1/Rp(j);
                    else
                        Rp(j) = inf;
                        indexJ(i,j+1) = i;
                    end

                end

                indexJ(i,1) = i;
                value(i,1) = sum(1 ./ Rp);

            end

        end

        function [value, sourceVector] = applyFixedTemperatureBC(obj, idx, value, sourceVector, resistances, TValue)

            % make idx as a row vector
            idx = idx(:)';
            Nidx = length(idx);
            NFN = obj.m.NFN;

            % calculate temperature value at facet centers -> tvfc
            if isa(TValue,'function_handle')
                tvfc = zeros(1,Nidx);
                for i = 1:Nidx
                    tvfc(i) = TValue(obj.m.facetCenter(idx(i),1),obj.m.facetCenter(idx(i),2),obj.m.facetCenter(idx(i),3));
                end
            else
                tvfc = TValue * ones(1,Nidx);
            end

            % loop over facets
            for i = 1:Nidx
                if obj.m.facets(idx(i),NFN+3)

                    eIndex = obj.m.facets(idx(i),NFN+3);
                    value(eIndex,1) = value(eIndex,1) + 1/resistances(eIndex,obj.m.facets(idx(i),NFN+4));
                    sourceVector(eIndex) = sourceVector(eIndex) + tvfc(i)/resistances(eIndex,obj.m.facets(idx(i),NFN+4));

                else

                    eIndex = obj.m.facets(idx(i),NFN+5);
                    value(eIndex,1) = value(eIndex,1) + 1/resistances(eIndex,obj.m.facets(idx(i),NFN+6));
                    sourceVector(eIndex) = sourceVector(eIndex) + tvfc(i)/resistances(eIndex,obj.m.facets(idx(i),NFN+6));

                end
            end

        end

        function [value, sourceVector] = applyConvectionBC(obj, idx, value, sourceVector, resistances, Tinf, HValue)

            % make idx as a row vector
            idx = idx(:)';
            Nidx = length(idx);
            um = obj.units.k_length;
            NFN = obj.m.NFN;

            % calculate h-coefficient value at facet centers -> hvfc
            if isa(HValue,'function_handle')
                hvfc = zeros(1,Nidx);
                for i = 1:Nidx
                    hvfc(i) = HValue(obj.m.facetCenter(idx(i),1),obj.m.facetCenter(idx(i),2),obj.m.facetCenter(idx(i),3));
                end
            else
                hvfc = HValue * ones(1,Nidx);
            end

            % loop over facets
            for i = 1:Nidx

                if obj.m.facets(idx(i),NFN+3)

                    eIndex = obj.m.facets(idx(i),NFN+3);
                    Gconv = 1/(resistances(eIndex,obj.m.facets(idx(i),NFN+4)) + 1/(hvfc(i) * obj.m.facetArea(idx(i)) * um^2));

                else

                    eIndex = obj.m.facets(idx(i),NFN+5);
                    Gconv = 1/(resistances(eIndex,obj.m.facets(idx(i),NFN+6)) + 1/(hvfc(i) * obj.m.facetArea(idx(i)) * um^2));

                end

                value(eIndex,1) = value(eIndex,1) + Gconv;
                sourceVector(eIndex) = sourceVector(eIndex)  + Tinf*Gconv;

            end

        end

        function sourceVector = applyHeatFluxBC(obj, idx, sourceVector, HFValue)

            % make idx as a row vector
            idx = idx(:)';
            Nidx = length(idx);

            % calculate heat flux at facet centers -> hffc
            if isa(HFValue,'function_handle')
                hffc = zeros(1,Nidx);
                for i = 1:Nidx
                    hffc(i) = HFValue(obj.m.facetCenter(idx(i),1),obj.m.facetCenter(idx(i),2),obj.m.facetCenter(idx(i),3));
                end
            else
                hffc = HFValue * ones(1,Nidx);
            end

            % check mesh type to find proper indices
            switch obj.meshType

                % tetrahedral mesh
                case 'thm'

                    for i = 1:Nidx

                        if obj.m.facets(idx(i),6)

                            eIndex = obj.m.facets(idx(i),6);
                            sourceVector(eIndex) = sourceVector(eIndex) + hffc(i)*obj.m.facetArea(idx(i));

                        else

                            eIndex = obj.m.facets(idx(i),8);
                            sourceVector(eIndex) = sourceVector(eIndex) + hffc(i)*obj.m.facetArea(idx(i));

                        end

                    end

                    % hexahedral mesh
                case {'pm','hhm'}

                    for i = 1:Nidx

                        if obj.m.facets(idx(i),7)

                            eIndex = obj.m.facets(idx(i),7);
                            sourceVector(eIndex) = sourceVector(eIndex) + hffc(i)*obj.m.facetArea(idx(i));

                        else

                            eIndex = obj.m.facets(idx(i),9);
                            sourceVector(eIndex) = sourceVector(eIndex) + hffc(i)*obj.m.facetArea(idx(i));

                        end

                    end

                otherwise
                    error('Mesh type is not supported.');
            end

        end

        function sourceVector = applyIHG(obj, mzName, sourceVector, IHGValue)

            idx = find(obj.m.ezi(:,obj.m.mzs.(mzName).zi));
            Nidx = length(idx);
            um = obj.units.k_length;

            % calculate internl heat generation at elment centers -> ihgec
            if isa(IHGValue,'function_handle')
                ihgec = zeros(1,Nidx);
                for i = 1:Nidx
                    ihgec(i) = IHGValue(obj.m.elementCenter(idx(i),1),obj.m.elementCenter(idx(i),2),obj.m.elementCenter(idx(i),3));
                end
            else
                ihgec = IHGValue * ones(1,Nidx);
            end

            for i = 1:Nidx

                sourceVector(idx(i)) = sourceVector(idx(i)) + ihgec(i) * obj.m.gev(idx(i)) * um^3;

            end

        end

        function gMatrix = applyContact(obj, idx, gMatrix, hValue)

            idx = idx(:)';
            um = obj.units.k_length;

            for index = idx

                % index of left element
                i = obj.m.facets(index,obj.m.NFN+3);

                % index of right element
                j = obj.m.facets(index,obj.m.NFN+5);

                % length of contact edge
                A = obj.m.facetArea(index) * um^2;

                if isa(hValue,'function_handle')
                    h = hValue(obj.m.facetCenter(index,1),obj.m.facetCenter(index,2),obj.m.facetCenter(index,3));
                else
                    h = hValue;
                end

                Rc = 1/(h*A);
                Rij_old = -1/gMatrix(i,obj.m.facets(index,obj.m.NFN+4) + 1);
                Rij_new = Rij_old + Rc;
                gMatrix(i,1) = gMatrix(i,1) - 1/Rij_old + 1/Rij_new;
                gMatrix(j,1) = gMatrix(j,1) - 1/Rij_old + 1/Rij_new;
                gMatrix(i,obj.m.facets(index,obj.m.NFN+4) + 1) = -1/Rij_new;
                gMatrix(j,obj.m.facets(index,obj.m.NFN+6) + 1) = -1/Rij_new;

            end

        end

    end

end
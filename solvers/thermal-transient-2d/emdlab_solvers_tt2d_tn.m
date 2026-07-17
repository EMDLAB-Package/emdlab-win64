% EMDLAB: Electrical Machines Design Laboratory
% a two-dimensional thermal-transient solver based on thermal netwrok

classdef emdlab_solvers_tt2d_tn < handle

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
        kxySV (:,:) double;

        % initial temperature assumption for iterative solver
        initialTemperature (1,1) double = 25;

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

        simTime (:,1); % simulation time vector
        simTimeStep (1,1) double = 1; % time step size for simulation
        simTimeStop (1,1) double = 60; % stop time for simulation

        % states
        isBeEvaluated (1,1) logical = false;
        isBnEvaluated (1,1) logical = false;
        isElementDataAssigned (1,1) logical = false;
        isResultsValid (1,1) logical = false;

        meshType (1,:) char = '';
        NR (1,1) double;

    end

    properties (Dependent = true)

        % number of excitations
        Nex (1,1) double;

    end

    methods
        %% Constructor and Destructor
        function obj = emdlab_solvers_tt2d_tn(m)

            % mesh pointer
            if isa(m, 'emdlab_m2d_qmdb')
                obj.meshType = 'qm';
                obj.NR = 4;
            elseif isa(m, 'emdlab_m2d_tmdb')
                obj.meshType = 'tm';
                obj.NR = 3;
                m.evalJIT;
            elseif isa(m, 'emdlab_m2d_qtmdb')
                obj.meshType = 'qtm';
            else
                error('Mesh class must be <emdlab_m2d_qmdb> or <emdlab_m2d_qmdb> or <emdlab_m2d_qtmdb>');
            end
            m.ggmesh;
            obj.m = m;

            % set default properties of mesh zones
            mzNames = obj.m.getMeshZoneNames;

            for mz = mzNames
                obj.setdp(mz);
            end

            % default values
            obj.depth = 1;
            %             obj.bcs = emdlab_bcs_scalarNodes('TL3');
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
        function setDepth(obj, value)
            obj.depth = value;
        end

        function y = getDepth(obj)
            y = obj.depth * obj.units.k_length;
        end
        %% solver properties for mesh zones

        function setdp(obj, mzName)
            % set default properties of a mesh zone

            obj.m.mzs.(mzName).props.isExcited = false;
            obj.m.mzs.(mzName).props.isCoilArm = false;
            obj.m.mzs.(mzName).props.isCoilMember = false;
            obj.m.mzs.(mzName).props.isMagnetized = false;
            obj.m.mzs.(mzName).props.initialTemperature = 25; 
            %             obj.makeFalse_isElementDataAssigned;

        end

        function setMeshZoneInitialTemperature(obj, mzName, value)
            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.m.mzs.(mzName).props.initialTemperature = value;
        end


        function addmz(obj, mzName, mzValue)

            obj.m.addmz(mzName, mzValue);
            obj.setdp(mzName);
            %             obj.makeFalse_isElementDataAssigned;

        end

        function removemz(obj, mzName)

            obj.m.removemz(mzName);
            %             obj.makeFalse_isElementDataAssigned;

        end

        function rotateMeshZone(obj, mzName, varargin)

            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.m.rotateMeshZone(mzName, varargin{:});

            if obj.m.mzs.(mzName).props.isMagnetized
                obj.m.mzs.(mzName).props.magnetization.rotate(varargin{:});
            end

        end

        %% Solver Functions
        function setSimulationStopTime(obj, value)
            obj.simTimeStop = value;
        end

        function setSimulationTimeStep(obj, value)
            obj.simTimeStep = value;
        end

        function setSolverIndex(obj, value)
            obj.solverIndex = value;
        end

        function assignElementData(obj)

            % check states
            if obj.isElementDataAssigned, return; end

            % allocation of memory
            obj.edata.ThermalConductivity = zeros(3, obj.m.Ne);
            obj.edata.MassDensity = zeros(1, obj.m.Ne);
            obj.edata.HeatCapacity = zeros(1, obj.m.Ne);

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

                obj.edata.MassDensity(1,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).MassDensity.value;
                obj.edata.HeatCapacity(1,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).HeatCapacity.value;

            end

            % change states
            obj.isElementDataAssigned = true;

        end

        function solve(obj)

            % assign element data
            obj.assignElementData;

            % calculate resistances
            resistances = zeros(obj.m.Ne,obj.NR);
            capacitances = zeros(obj.m.Ne,1);
            z = obj.getDepth;

            % loop over mesh zones
            mzNames = obj.m.getMeshZoneNames;
            ux = zeros(obj.m.Ne,3);
            uy = zeros(obj.m.Ne,3);            
            for mz = mzNames
                mzptr = obj.m.mzs.(mz);                
                ux(obj.m.ezi(:,mzptr.zi),:) = repmat(obj.m.cs.(mzptr.orientation).ux,mzptr.Ne,1);
                uy(obj.m.ezi(:,mzptr.zi),:) = repmat(obj.m.cs.(mzptr.orientation).uy,mzptr.Ne,1);
            end

            % store source term due to anisotropy
            obj.kxySV = zeros(obj.m.Ne,obj.NR);

            % element centers
            elementCenter = obj.m.getCenterOfElements;
            edgeCenter = obj.m.getCenterOfEdges;
            elm = zeros(obj.m.Ne,obj.NR);

            % loop over elements to calculate conductances of each element
            for i = 1:obj.m.Ne

                kx = obj.edata.ThermalConductivity(1,i);
                ky = obj.edata.ThermalConductivity(2,i);

                for j = 1:obj.NR

                    % direction vector from cell center i to cell center j
                    if obj.m.nbs(i,j)
                        uij = elementCenter(obj.m.nbs(i,j), :) - elementCenter(i,:);
                        die = 0.5 * norm(uij);
                    else
                        uij = edgeCenter(abs(obj.m.elements(i,j)),:) - elementCenter(i,:);
                        die = norm(uij);
                    end                    
                    uij = uij / norm(uij);

                    kij = kx * (ux(i,1:2) * uij').^2 + ky * (uy(i,1:2) * uij').^2;
                    obj.kxySV(i,j) = (kx-ky) * (ux(i,1:2) * uij') * (uy(i,1:2) * uij');

                    elm(i,j) = abs(uij * obj.m.nEdges(abs(obj.m.elements(i,j)),:)');

                    % calculation of conductances
                    resistances(i,j) = die/(kij*z*obj.m.el(i,j) * elm(i,j));

                    

                end

                capacitances(i) = obj.edata.MassDensity(i) * obj.edata.HeatCapacity(i) * obj.m.gea(i) * z * obj.units.k_length^2;

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

            obj.simTime = linspace(0, obj.simTimeStop, ceil(obj.simTimeStop/obj.simTimeStep))';
            Nt = length(obj.simTime);

            for mz = mzNames
                obj.results.(mz).T = zeros(Nt,3);
            end

            C_Matrix = sparse(1:obj.m.Ne, 1:obj.m.Ne, capacitances)/obj.simTime(2);
            G_maxtrix = G_maxtrix + C_Matrix;

            % set initial temperature
            T_old = zeros(obj.m.Ne,1);
            for mz = mzNames
                T_old(obj.m.ezi(:,obj.m.mzs.(mz).zi)) = obj.m.mzs.(mz).props.initialTemperature;
                obj.results.(mz).T(1,:) = obj.m.mzs.(mz).props.initialTemperature;
            end
            
            obj.results.T = T_old;

            % loop over time
            for k = 2:Nt
                
                obj.results.T = G_maxtrix\(sourceVector + C_Matrix * T_old);
                obj.evalTn;

                obj.solverHistory.relativeError = [];

                if ~obj.edata.areAllTemperatureIndependent || ~obj.edata.areAllIsotropic || ~obj.edata.areAllHomogeneous

                    releativeError = 1e-3;
                    maxIteration = 200;
                    iter = 0;
                    err = inf;

                    % iterative loop for solver
                    while iter<maxIteration && err>releativeError

                        sourceVectorU = (sourceVector + C_Matrix * T_old);

                        % handle anisotropy
                        if ~obj.edata.areAllIsotropic
                            for i = 1:obj.m.Ne
                                Tij = obj.results.Tn(obj.m.cl(i,[1:obj.NR,1]));
                                for j = 1:obj.NR
                                    if ~obj.m.bedges(abs(obj.m.elements(i,j)))
                                        sourceVectorU(i) = sourceVectorU(i) + z * ...
                                            obj.kxySV(i,j) * (Tij(j) - Tij(j+1)) * elm(i,j);
                                    end
                                end
                            end
                        end

                        Told = obj.results.T;

                        % solver iteration
                        obj.results.T = G_maxtrix\sourceVectorU;
                        obj.evalTn;

                        iter = iter + 1;
                        err = norm(obj.results.T-Told,2)/norm(obj.results.T,2);
                        obj.solverHistory.relativeError(end+1) = err;

                        fprintf('Iteration #%03d, Relative Error = %.2e\n', iter, err);

                    end

                end

                T_old = obj.results.T;

                for mz = mzNames
                    obj.results.(mz).T(k,:) = [min(T_old), max(T_old), mean(T_old)];
                end

                fprintf('simulation time step %.2e completed.\n', obj.simTime(k));

            end

            % calculate cross heat at edges
            Nedges = size(obj.m.edges,1);
            obj.results.qe = zeros(Nedges,1);
            for k = 1:Nedges
                if ~obj.m.bedges(k)
                    % index of left element
                    i = obj.m.edges(k,5);
                    % index of right element
                    j = obj.m.edges(k,7);
                    % net heat passing k-th internal edge
                    obj.results.qe(k) = (obj.results.T(j) - obj.results.T(i)) * G_maxtrix(i,j);
                end
            end

            for k = 1:numel(exNames)

                exptr = obj.excitations.(exNames(k));
                switch exptr.type
                    
                    case 'fixed-temperature'
                        if isa(exptr.value, 'function_handle')
                                Tb = exptr.value(edgeCenter(exptr.idx,1),edgeCenter(exptr.idx,2));
                            else
                                Tb = exptr.value*ones(length(exptr.idx),1);
                        end

                        j = 0;
                        for i = exptr.idx(:)'
                            j = j+1;

                            if obj.m.edges(i,4)
                                obj.results.qe(i) = (Tb(j) - obj.results.T(obj.m.edges(i,7))) / ...
                                    resistances(obj.m.edges(i,7), obj.m.edges(i,8));
                            else
                                obj.results.qe(i) = (obj.results.T(obj.m.edges(i,5)) - Tb(j)) / ...
                                    resistances(obj.m.edges(i,5), obj.m.edges(i,6));
                            end
                        end

                    case 'convection'
                        if isa(exptr.hValue, 'function_handle')
                                hValue = exptr.hValue(edgeCenter(exptr.idx,1),edgeCenter(exptr.idx,2));
                            else
                                hValue = exptr.hValue*ones(length(exptr.idx),1);
                        end

                        j = 0;
                        for i = exptr.idx(:)'
                            j = j+1;

                            if obj.m.edges(i,4)
                                obj.results.qe(i) = (exptr.Tinf - obj.results.T(obj.m.edges(i,7))) / ...
                                    (resistances(obj.m.edges(i,7), obj.m.edges(i,8)) + 1/(hValue(j)*z*obj.m.edgeLength(i)*obj.units.k_length));
                            else
                                obj.results.qe(i) = (obj.results.T(obj.m.edges(i,5)) - exptr.Tinf) / ...
                                    (resistances(obj.m.edges(i,5), obj.m.edges(i,6)) + 1/(hValue(j)*z*obj.m.edgeLength(i)*obj.units.k_length));
                            end
                        end

                    case 'radiation'
                        [value, sourceVector] = applyConvectionBC(obj, exptr.idx, value, sourceVector, resistances, exptr.Tinf, exptr.hValue);
                    
                    case 'heat-flux'
                        if isa(exptr.value, 'function_handle')
                            qValue = exptr.value(edgeCenter(exptr.idx,1),edgeCenter(exptr.idx,2));
                        else
                            qValue = exptr.value*ones(length(exptr.idx),1);
                        end

                        j = 0;
                        for i = exptr.idx(:)'
                            j = j+1;

                            if obj.m.edges(i,4)
                                obj.results.qe(i) = qValue(j)*z*obj.m.edgeLength(i)*obj.units.k_length;
                            else
                                obj.results.qe(i) = -qValue(j)*z*obj.m.edgeLength(i)*obj.units.k_length;
                            end
                        end

                end

            end

        end

        function y = calculateNetHeatCrossingBoundaryEdges(obj, idx)
            y = 0;
            for i = idx(:)'
                if obj.m.edges(i,4)
                    y = y + obj.results.qe(i);
                else
                    y = y - obj.results.qe(i);
                end
            end
        end

        function y = calculateNetHeatCrossingContact(obj, cName)
            y = 0;
            for i = idx(:)'
                if obj.m.edges(i,4)
                    y = y + obj.results.qe(i);
                else
                    y = y - obj.results.qe(i);
                end
            end
        end

        function evalTn(obj)

            obj.results.Tn = zeros(obj.m.Nn,1);
            ec = obj.m.getCenterOfElements;
            denum = zeros(obj.m.Nn,1);

            % loop over elements
            for i = 1:obj.m.Ne
                for j = 1:obj.NR
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

            gea_tmp = repmat(obj.m.gea,obj.NR,1);
            index = obj.m.cl';
            obj.results.Tsmooth = full(diag(sparse(index, index, repmat(obj.results.T',obj.NR,1) .* gea_tmp)) ./ ...
            diag(sparse(index, index, gea_tmp)));
            
            % apply excitation conditions
            exNames = obj.getExcitationNames;
            for i = 1:numel(exNames)
                exptr = obj.excitations.(exNames(i));
                if strcmpi(exptr.type,'fixed-temperature')
                    for j = reshape(exptr.idx, 1, [])
                        if isa(exptr.value, 'function_handle')
                            pts = obj.m.nodes(obj.m.edges(j,1:2),:);
                            obj.results.Tsmooth(obj.m.edges(j,1)) = exptr.value(pts(1,1),pts(1,2));
                            obj.results.Tsmooth(obj.m.edges(j,2)) = exptr.value(pts(2,1),pts(2,2));
                        else
                            obj.results.Tsmooth(obj.m.edges(j,1:2)) = exptr.value;
                        end
                    end                     
                end
            end            

        end

        %% Excitation Definitions
        function addFixedTemperatureBC(obj, exName, idx, value)

            if isempty(idx), error('Empty edge indices.'); end

            exName = obj.checkExcitationNonExistence(exName);
            obj.excitations.(exName).type = 'fixed-temperature';
            obj.excitations.(exName).idx = idx;
            obj.excitations.(exName).value = value;

        end

        function addConvectionBC(obj, exName, idx, hValue, Tinf)

            if isempty(idx), error('Empty edge indices.'); end

            exName = obj.checkExcitationNonExistence(exName);
            obj.excitations.(exName).type = 'convection';
            obj.excitations.(exName).idx = idx;
            obj.excitations.(exName).hValue = hValue;
            obj.excitations.(exName).Tinf = Tinf;

        end

        function addRadiationBC(obj, exName, idx, eValue, Tsurr)

            if isempty(idx), error('Empty edge indices.'); end

            exName = obj.checkExcitationNonExistence(exName);
            obj.excitations.(exName).type = 'radiation';
            obj.excitations.(exName).idx = idx;
            obj.excitations.(exName).eValue = eValue;
            obj.excitations.(exName).Tsurr = Tsurr;

        end

        function addHeatFluxBC(obj, exName, idx, value)

            if isempty(idx), error('Empty edge indices.'); end

            exName = obj.checkExcitationNonExistence(exName);
            obj.excitations.(exName).type = 'heat-flux';
            obj.excitations.(exName).idx = idx;
            obj.excitations.(exName).value = value;

        end

        function addHeatFlowBC(obj, exName, idx, value)

            if isempty(idx), error('Empty edge indices.'); end

            exName = obj.checkExcitationNonExistence(exName);
            obj.excitations.(exName).type = 'heat-flux';
            obj.excitations.(exName).idx = idx;
            area = sum(obj.m.edgeLength(idx)) * obj.getDepth * obj.units.k_length;
            obj.excitations.(exName).value = value/area;

        end

        function addInternalHeatSource(obj, exName, mzName, value, unit)

            exName = obj.checkExcitationNonExistence(exName);
            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.excitations.(exName).type = 'internal-heat-source';
            obj.excitations.(exName).mzName = mzName;
            
            z = obj.getDepth;
            um = obj.units.k_length;
            Area = sum(obj.m.gea(obj.m.ezi(:,obj.m.mzs.(mzName).zi))) * um^2;

            % set default unit
            if nargin < 5, unit = 'W'; end

            unit = lower(erase(unit, ' '));
            switch unit

                case 'w'
                    value = value/(z*Area);

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
            obj.contacts.(cName).idx = obj.m.getEdgeIndicesOnContact(mz1Name,mz2Name);
            obj.contacts.(cName).mz1Name = mz1Name;
            obj.contacts.(cName).mz2Name = mz2Name;
            obj.contacts.(cName).hValue = hValue;

        end

        function area = getContactArea(obj, mz1Name, mz2Name)
            idx = obj.m.getEdgeIndicesOnContact(mz1Name,mz2Name);
            area = obj.getDepth * sum(obj.m.edgeLength(idx)) * obj.units.k_length;
        end

        %% Visualization Functions
        function varargout = plotAverageTemperature(obj, N, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            if nargin<2, N=10; end

            patch('faces', obj.m.cl, 'Vertices', obj.m.nodes, 'FaceVertexCData',obj.results.T, ...
                'FaceColor','flat', 'edgecolor', 'none');
            colormap(jet(N));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Average Temperature [C]';
            climits = clim; 
            cb.Ticks = fix(linspace(climits(1), climits(2), 10)*100)/100;
            cb.Ticks(1) = cb.Ticks(1) + 0.01;

            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotTemperature(obj, N, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            if nargin<2, N=10; end
            obj.evalTSmooth;

            patch('faces', obj.m.cl, 'Vertices', obj.m.nodes, 'FaceVertexCData', obj.results.Tsmooth, ...
                'FaceColor','interp', 'edgecolor', 'none');
            colormap(jet(N));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Temperature [C]';
            climits = clim; 
            cb.Ticks = fix(linspace(climits(1), climits(2), 10)*100)/100;
            cb.Ticks(1) = cb.Ticks(1) + 0.01;

            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotTemperatureTn(obj, N, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            if nargin<2, N=10; end
            obj.evalTn;

            patch('faces', obj.m.cl, 'Vertices', obj.m.nodes, 'FaceVertexCData', obj.results.Tn, ...
                'FaceColor','interp', 'edgecolor', 'none');
            colormap(jet(N));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Temperature [C]';
            climits = clim; 
            cb.Ticks = fix(linspace(climits(1), climits(2), 10)*100)/100;
            cb.Ticks(1) = cb.Ticks(1) + 0.01;

            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function y = getAverageTemperature(obj)
            y = obj.m.gea * obj.results.T / sum(obj.m.gea);
        end

        function showContact(obj, varargin)

            if numel(varargin) == 1
                cName = obj.checkContactExistence(varargin{1});
                idx = obj.contacts.(cName).idx;
            elseif numel(varargin) == 2
                mz1Name = obj.m.checkMeshZoneExistence(varargin{1});
                mz2Name = obj.m.checkMeshZoneExistence(varargin{2});
                idx = obj.m.getEdgeIndicesOnContact(mz1Name, mz2Name);
            else
                error('Wrong inputs.');
            end

            obj.m.showEdges(idx);

        end

        function plotMeshZoneTemperatureVsTime(obj, mzName)

            figure; box on; hold on;
             plot(obj.simTime, obj.results.(mzName).T, 'LineWidth',1.5);
             xlabel('Time [s]')
             ylabel('Temperature [C]');
             legend('Min', 'Max', 'Average');

        end

        % show center to edge connections
        function varargout = plotThermalNetwork(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});

            mzNames = obj.m.getMeshZoneNames;
            ecolor = zeros(obj.m.Ne, 3);

            for i = 1:obj.m.Nmzs
                mzptr = obj.m.mzs.(mzNames{i});
                ecolor(obj.m.ezi(:, i), 1) = mzptr.color(1);
                ecolor(obj.m.ezi(:, i), 2) = mzptr.color(2);
                ecolor(obj.m.ezi(:, i), 3) = mzptr.color(3);
            end

            pts = [obj.m.getCenterOfElements; obj.m.getCenterOfEdges];
            cl1_tmp = repmat(1:obj.m.Ne,obj.NR,1);
            cl2_tmp = abs(obj.m.elements(:,1:obj.NR))' + obj.m.Ne;
            cl3_tmp = [cl1_tmp(:),cl2_tmp(:)];

            patch('Faces', obj.m.cl(:, 1:obj.NR), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'flat', 'FaceVertexCData', ecolor, ...
                'FaceAlpha', 1, 'EdgeColor', 'k', 'linewidth', 1, 'parent', ax);

            P = pts;

            cl = cl3_tmp;

            Xall = [];
            Yall = [];

            for e = 1:size(cl,1)
                i = cl(e,1);
                j = cl(e,2);

                p1 = P(i,:);
                p2 = P(j,:);

                d = p2 - p1;
                L = norm(d);
                if L == 0
                    continue
                end

                ex = d / L;
                ey = [-ex(2), ex(1)];

                leadFrac = 0.15;
                leadLen  = leadFrac * L;
                amp      = 0.08 * L;
                nzig     = 6;

                xzig = linspace(leadLen, L-leadLen, 2*nzig+1);
                yzig = zeros(size(xzig));
                yzig(2:2:end-1) = amp;
                yzig(3:2:end-1) = -amp;

                xloc = [0, leadLen, xzig, L];
                yloc = [0, 0,       yzig, 0];

                pts = p1 + xloc.'*ex + yloc.'*ey;

                Xall = [Xall; pts(:,1); NaN];
                Yall = [Yall; pts(:,2); NaN];
            end

            
            patch(Xall, Yall, 'w', ...
                'EdgeColor', 'b', ...
                'FaceColor', 'none', ...
                'LineWidth', 1);

%             plot(P(:,1), P(:,2), 'ro', 'MarkerFaceColor', 'r');
%             text(P(:,1), P(:,2), compose(' %d',1:size(P,1)));


%             patch('Faces', cl3_tmp, 'Vertices', pts, ...
%                 'FaceColor', 'none', 'EdgeColor', 'b', 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

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
                for j = 1:obj.NR

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

            if iscolumn(idx), idx = idx'; end
            edgeCenters = obj.m.getCenterOfEdges;

            for index = idx

                if isa(TValue,'function_handle')
                    tvbc = TValue(edgeCenters(index,1),edgeCenters(index,2));
                else
                    tvbc = TValue;
                end

                if obj.m.edges(index,5)

                    eIndex = obj.m.edges(index,5);
                    value(eIndex,1) = value(eIndex,1) + 1/resistances(eIndex,obj.m.edges(index,6));
                    sourceVector(eIndex) = sourceVector(eIndex) + tvbc/resistances(eIndex,obj.m.edges(index,6));

                else

                    eIndex = obj.m.edges(index,7);
                    value(eIndex,1) = value(eIndex,1) + 1/resistances(eIndex,obj.m.edges(index,8));
                    sourceVector(eIndex) = sourceVector(eIndex) + tvbc/resistances(eIndex,obj.m.edges(index,8));

                end

            end

        end

        function [value, sourceVector] = applyConvectionBC(obj, idx, value, sourceVector, resistances, Tinf, HValue)

            if iscolumn(idx), idx = idx'; end
            edgeCenters = obj.m.getCenterOfEdges;
            z = obj.getDepth;
            um = obj.units.k_length;

            for index = idx

                if isa(HValue,'function_handle')
                    hvbc = HValue(edgeCenters(index,1),edgeCenters(index,2));
                else
                    hvbc = HValue;
                end

                if obj.m.edges(index,5)

                    eIndex = obj.m.edges(index,5);
                    if ismember(obj.m.edges(index,6),[1,3])
                        Gconv = 1/(resistances(eIndex,1)/2 + 1/(hvbc*z*um*obj.m.el(eIndex,obj.m.edges(index,6))));
                    else
                        Gconv = 1/(resistances(eIndex,2)/2 + 1/(hvbc*z*um*obj.m.el(eIndex,obj.m.edges(index,6))));
                    end

                else

                    eIndex = obj.m.edges(index,7);
                    if ismember(obj.m.edges(index,8),[1,3])
                        Gconv = 1/(resistances(eIndex,1)/2 + 1/(hvbc*z*um*obj.m.el(eIndex,obj.m.edges(index,8))));
                    else
                        Gconv = 1/(resistances(eIndex,2)/2 + 1/(hvbc*z*um*obj.m.el(eIndex,obj.m.edges(index,8))));
                    end

                end

                value(eIndex,1) = value(eIndex,1) + Gconv;
                sourceVector(eIndex) = sourceVector(eIndex)  + Tinf*Gconv;

            end

        end

        function sourceVector = applyHeatFluxBC(obj, idx, sourceVector, HFValue)

            if iscolumn(idx), idx = idx'; end
            edgeCenters = obj.m.getCenterOfEdges;
            z = obj.getDepth;
            um = obj.units.k_length;

            for index = idx

                if isa(HFValue,'function_handle')
                    hfbc = HFValue(edgeCenters(index,1),edgeCenters(index,2));
                else
                    hfbc = HFValue;
                end

                if obj.m.edges(index,5)

                    eIndex = obj.m.edges(index,5);
                    sourceVector(eIndex) = sourceVector(eIndex) + hfbc*obj.m.el(eIndex,obj.m.edges(index,6))*um*z;

                else

                    eIndex = obj.m.edges(index,7);
                    sourceVector(eIndex) = sourceVector(eIndex) + hfbc*obj.m.el(eIndex,obj.m.edges(index,8))*um*z;

                end

            end

        end

        function sourceVector = applyIHG(obj, mzName, sourceVector, IHGValue)

            idx = 1:obj.m.Ne;
            idx = idx(obj.m.ezi(:,obj.m.mzs.(mzName).zi));
            z = obj.getDepth;
            um = obj.units.k_length;

            for index = idx

                if isa(IHGValue,'function_handle')
                    ihgbc = IHGValue(edgeCenters(index,1),edgeCenters(index,2));
                else
                    ihgbc = IHGValue;
                end

                sourceVector(index) = sourceVector(index) + ihgbc * obj.m.gea(index) * z * um^2;

            end

        end

        function gMatrix = applyContact(obj, idx, gMatrix, hValue)

            if iscolumn(idx), idx = idx'; end
            edgeCenters = obj.m.getCenterOfEdges;
            z = obj.getDepth;
            um = obj.units.k_length;

            for index = idx

                % index of left element
                i = obj.m.edges(index,5);

                % index of right element
                j = obj.m.edges(index,7);

                % length of contact edge
                el = obj.m.el(i,obj.m.edges(index,6));

                % area of contact edge
                A = el*z*um;

                if isa(hValue,'function_handle')
                    h = hValue(edgeCenters(index,1),edgeCenters(index,2));
                else
                    h = hValue;
                end

                Rc = 1/(h*A);
                Rij_old = -1/gMatrix(i,obj.m.edges(index,6) + 1);
                Rij_new = Rij_old + Rc;
                gMatrix(i,1) = gMatrix(i,1) - 1/Rij_old + 1/Rij_new;
                gMatrix(j,1) = gMatrix(j,1) - 1/Rij_old + 1/Rij_new;
                gMatrix(i,obj.m.edges(index,6) + 1) = -1/Rij_new;
                gMatrix(j,obj.m.edges(index,8) + 1) = -1/Rij_new;

            end

        end

    end

end
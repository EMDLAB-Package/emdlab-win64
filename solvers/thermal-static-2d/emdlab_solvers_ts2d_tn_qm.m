% EMDLAB: Electrical Machines Design Laboratory
% a two-dimensional thermal-static solver based on TN for QM

classdef emdlab_solvers_ts2d_tn_qm < handle

    properties (SetAccess = protected)

        % solver mesh
        m (1,1) emdlab_m2d_qmdb;

        % internal & boundary conditions
        excitations (1,1) struct;

        % contacts
        contacts (1,1) struct;

        % elements data
        edata (1,1) struct;

        % results
        results (1,1) struct;

        kxySV (:,1) double;

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

        % number of excitations
        Nex (1,1) double;

        % number of coils
        Ncoils (1,1) double;

    end

    methods
        %% Constructor and Destructor
        function obj = emdlab_solvers_ts2d_tn_qm(m)

            % mesh pointer
            if isa(m, 'emdlab_m3d_qmdb')
                error('Mesh class must be <emdlab_m2d_qmdb');
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
            %             obj.makeFalse_isElementDataAssigned;

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
        function setSolverIndex(obj, value)
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
            resistances = zeros(obj.m.Ne,2);
            z = obj.getDepth;

            ec = obj.m.getCenterOfEdges;
            xec = ec(:,1);
            yec = ec(:,2);
            xec = xec(abs(obj.m.elements(:,1:4)));
            yec = yec(abs(obj.m.elements(:,1:4)));

            ux = zeros(obj.m.Ne,3);
            uy = zeros(obj.m.Ne,3);

            mzNames = obj.m.getMeshZoneNames;
            for mz = mzNames
                mzptr = obj.m.mzs.(mz);
                
                ux(obj.m.ezi(:,mzptr.zi),:) = repmat(obj.m.cs.(mzptr.orientation).ux,mzptr.Ne,1);
                uy(obj.m.ezi(:,mzptr.zi),:) = repmat(obj.m.cs.(mzptr.orientation).uy,mzptr.Ne,1);
            end

            obj.kxySV = zeros(obj.m.Ne,1);
            % loop over elements to calculate conductances of each element
            for i = 1:obj.m.Ne

                kx = obj.edata.ThermalConductivity(1,i);
                ky = obj.edata.ThermalConductivity(2,i);

                u13 = [xec(i,3),yec(i,3),0] - [xec(i,1),yec(i,1),0];
                u24 = [xec(i,4),yec(i,4),0] - [xec(i,2),yec(i,2),0];
                u13 = u13 / norm(u13);
                u24 = u24 / norm(u24);

                k13 = kx * dot(ux(i,:),u13).^2 + ky * dot(uy(i,:),u13).^2;
                k24 = kx * dot(ux(i,:),u24).^2 + ky * dot(uy(i,:),u24).^2;
                obj.kxySV(i) = (kx-ky) * (dot(ux(i,:),u13)*dot(uy(i,:),u13) - dot(ux(i,:),u24)*dot(uy(i,:),u24))/2;


                % calculation of conductances
                resistances(i,1) = (obj.m.el(i,2) + obj.m.el(i,4))/(k13*z*(obj.m.el(i,1) + obj.m.el(i,3)));
                resistances(i,2) = (obj.m.el(i,1) + obj.m.el(i,3))/(k24*z*(obj.m.el(i,2) + obj.m.el(i,4)));

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
            obj.evalTSmooth;
            sourceVectorU = sourceVector;

            if ~obj.edata.areAllTemperatureIndependent || ~obj.edata.areAllIsotropic || ~obj.edata.areAllHomogeneous

                releativeError = 1e-3;
                maxIteration = 50;
                iter = 0;
                err = inf;
                
                % iterative loop for solver
                while iter<maxIteration && err>releativeError

                    % handle anisotropy
                    if ~obj.edata.areAllIsotropic
                        for i = 1:obj.m.Ne
                            T1 = obj.results.Tsmooth(obj.m.cl(i,1));
                            T2 = obj.results.Tsmooth(obj.m.cl(i,2));
                            T3 = obj.results.Tsmooth(obj.m.cl(i,3));
                            T4 = obj.results.Tsmooth(obj.m.cl(i,4));
                            sourceVectorU(i) = sourceVector(i) + obj.kxySV(i)*z*  ...
                                (-(T2-T1)+(T3-T2)-(T4-T3)+(T1-T4));
                        end
                    end

                    Told = obj.results.T;
                    % solver iteration
                    obj.results.T = G_maxtrix\sourceVectorU;
                    obj.evalTSmooth;

                    iter = iter + 1;
                    err = norm(obj.results.T-Told,2)/norm(obj.results.T,2);

                end

            end

        end

        function evalTSmooth(obj)

            gea_tmp = repmat(obj.m.gea,4,1);
            index = obj.m.cl';
            obj.results.Tsmooth = full(diag(sparse(index, index, repmat(obj.results.T',4,1) .* gea_tmp)) ./ ...
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

        function y = getAverageTemperature(obj)
            y = obj.m.gea * obj.results.T / sum(obj.m.gea);
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
            cl1_tmp = repmat(1:obj.m.Ne,4,1);
            cl2_tmp = abs(obj.m.elements(:,1:4))' + obj.m.Ne;
            cl3_tmp = [cl1_tmp(:),cl2_tmp(:)];

            patch('Faces', obj.m.cl(:, [1,2,3,4]), 'Vertices', obj.m.nodes, ...
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

            indexI = zeros(obj.m.Ne,5);
            indexJ = zeros(obj.m.Ne,5);
            value = zeros(obj.m.Ne,5);
            sourceVector = zeros(obj.m.Ne,1);

            % loop over elements to construct G_matrix
            for i = 1:obj.m.Ne

                % constructing i-th row
                indexI(i,:) = i;

                % conductance between cell P and 1 nb
                if obj.m.nbs(i,1)
                    if ismember(abs(obj.m.elements(i,1)),abs(obj.m.elements(obj.m.nbs(i,1),[1,3])))
                        Rp1 = resistances(i,1)/2 + resistances(obj.m.nbs(i,1),1)/2;
                    else
                        Rp1 = resistances(i,1)/2 + resistances(obj.m.nbs(i,1),2)/2;
                    end
                    indexJ(i,2) = obj.m.nbs(i,1);
                    value(i,2) = -1/Rp1;
                else
                    Rp1 = inf;
                    indexJ(i,2) = i;
                end

                % conductance between cell P and 2 nb
                if obj.m.nbs(i,2)
                    if ismember(abs(obj.m.elements(i,2)),abs(obj.m.elements(obj.m.nbs(i,2),[1,3])))
                        Rp2 = resistances(i,2)/2 + resistances(obj.m.nbs(i,2),1)/2;
                    else
                        Rp2 = resistances(i,2)/2 + resistances(obj.m.nbs(i,2),2)/2;
                    end
                    indexJ(i,3) = obj.m.nbs(i,2);
                    value(i,3) = -1/Rp2;
                else
                    Rp2 = inf;
                    indexJ(i,3) = i;
                end

                % conductance between cell P and 3 nb
                if obj.m.nbs(i,3)
                    if ismember(abs(obj.m.elements(i,3)),abs(obj.m.elements(obj.m.nbs(i,3),[1,3])))
                        Rp3 = resistances(i,1)/2 + resistances(obj.m.nbs(i,3),1)/2;
                    else
                        Rp3 = resistances(i,1)/2 + resistances(obj.m.nbs(i,3),2)/2;
                    end
                    indexJ(i,4) = obj.m.nbs(i,3);
                    value(i,4) = -1/Rp3;
                else
                    Rp3 = inf;
                    indexJ(i,4) = i;
                end

                % conductance between cell P and 4 nb
                if obj.m.nbs(i,4)
                    if ismember(abs(obj.m.elements(i,4)),abs(obj.m.elements(obj.m.nbs(i,4),[1,3])))
                        Rp4 = resistances(i,2)/2 + resistances(obj.m.nbs(i,4),1)/2;
                    else
                        Rp4 = resistances(i,2)/2 + resistances(obj.m.nbs(i,4),2)/2;
                    end
                    indexJ(i,5) = obj.m.nbs(i,4);
                    value(i,5) = -1/Rp4;
                else
                    Rp4 = inf;
                    indexJ(i,5) = i;
                end

                indexJ(i,1) = i;
                value(i,1) = 1/Rp1 + 1/Rp2 + 1/Rp3 + 1/Rp4;

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
                    if ismember(obj.m.edges(index,6),[1,3])
                        value(eIndex,1) = value(eIndex,1) + 2/resistances(eIndex,1);
                        sourceVector(eIndex) = sourceVector(eIndex) + tvbc*2/resistances(eIndex,1);
                    else
                        value(eIndex,1) = value(eIndex,1) + 2/resistances(eIndex,2);
                        sourceVector(eIndex) = sourceVector(eIndex) + tvbc*2/resistances(eIndex,2);
                    end

                else

                    eIndex = obj.m.edges(index,7);
                    if ismember(obj.m.edges(index,8),[1,3])
                        value(eIndex,1) = value(eIndex,1) + 2/resistances(eIndex,1);
                        sourceVector(eIndex) = sourceVector(eIndex) + tvbc*2/resistances(eIndex,1);
                    else
                        value(eIndex,1) = value(eIndex,1) + 2/resistances(eIndex,2);
                        sourceVector(eIndex) = sourceVector(eIndex) + tvbc*2/resistances(eIndex,2);
                    end

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
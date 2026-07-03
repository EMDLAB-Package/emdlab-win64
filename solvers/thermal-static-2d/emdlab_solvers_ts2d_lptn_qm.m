% EMDLAB: Electrical Machines Design Laboratory
% a two-dimensional thermal-static solver based on LPTN for QM

classdef emdlab_solvers_ts2d_lptn_qm < handle

    properties (SetAccess = protected)

        % solver mesh
        m (1,1) emdlab_m2d_qmdb;

        % internal & boundary conditions
        excitations (1,1) struct;

        % elements data
        edata (1,1) struct;

        % results
        results (1,1) struct;

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
        function obj = emdlab_solvers_ts2d_lptn_qm(m)

            % mesh pointer
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
            obj.edata.areAllLinear = true;

            % loop over mesh zones
            for i = 1:obj.m.Nmzs

                % get pointer to mesh zone
                mzptr = obj.m.mzs.(mzNames{i});

                % assigning thermal conductivities
                if obj.m.mts.(mzptr.material).ThermalConductivity.isIsotropic

                    obj.edata.ThermalConductivity(:,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ThermalConductivity.value;

                else

                    obj.edata.ThermalConductivity(:,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ThermalConductivity.value;

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

            % loop over elements to calculate conductances of each element
            for i = 1:obj.m.Ne

                k = obj.edata.ThermalConductivity(1,i);

                % calculation of conductances
                resistances(i,1) = (obj.m.el(i,2) + obj.m.el(i,4))/(k*z*(obj.m.el(i,1) + obj.m.el(i,3)));
                resistances(i,2) = (obj.m.el(i,1) + obj.m.el(i,3))/(k*z*(obj.m.el(i,2) + obj.m.el(i,4)));

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

            % construct global matrices and solve
            G_maxtrix = sparse(indexI, indexJ, value);
            obj.results.T = G_maxtrix\sourceVector;

        end

        function evalTSmooth(obj)

            obj.results.Tsmooth = zeros(obj.m.Nn,1);
            areaAttachedToNode = zeros(obj.m.Nn,1);
            for i = 1:obj.m.Ne
                for j = 1:4
                    obj.results.Tsmooth(obj.m.cl(i,j)) = obj.results.Tsmooth(obj.m.cl(i,j)) + ...
                        obj.results.T(i) * obj.m.gea(i);
                    areaAttachedToNode(obj.m.cl(i,j)) = areaAttachedToNode(obj.m.cl(i,j)) + obj.m.gea(i);
                end
            end
            obj.results.Tsmooth = obj.results.Tsmooth ./ areaAttachedToNode;
            
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

        function addInternalHeatSource(obj, exName, mzName, value)

            exName = obj.checkExcitationNonExistence(exName);
            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.excitations.(exName).type = 'internal-heat-source';
            obj.excitations.(exName).mzName = mzName;
            obj.excitations.(exName).value = value;

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

    end

end
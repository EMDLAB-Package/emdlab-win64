% EMDLAB: Electrical Machines Design Laboratory
% Prism mesh data base (3D element)

classdef emdlab_m3d_pmdb < handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_m3d_xmdb

    properties (Dependent = true)

        % flags for element type
        isPRISM6 (1,1) logical;

    end

    methods
        %% Constructor and Destructor
        function obj = emdlab_m3d_pmdb(varargin)
            % add default material
            obj.addMaterial('air');
        end

        function addMeshZone(obj, varargin)
            % adding a new mesh zone to mesh database

            % ipnut arguments check
            switch nargin

                case 2
                    % you can pass an <emdlab_m3d_pmz> class without a name
                    % a default name will be considered in this case
                    if ~isa(varargin{1}, 'emdlab_m3d_pmz')
                        error('Mesh zone class must be <emdlab_m3d_pmz>.');
                    end
                    mzName = obj.getDefaultMeshZoneName;
                    mzptr = varargin{1};

                    
                case 3
                    % you can pass an <emdlab_m3d_pmz> mesh zone with a specified name
                    if ischar(varargin{1}) || isstring(varargin{1})

                        mzName = obj.checkMeshZoneNonExistence(varargin{1});
                        if ~isa(varargin{2}, 'emdlab_m3d_pmz')
                            error('Mesh zone class must be <emdlab_m3d_hhmz>.');
                        end
                        mzptr = varargin{2};

                    elseif isnumeric(varargin{1}) && isnumeric(varargin{2})
                        % add a new mesh zone by direct passing connectivity list and points
                        % considering a default mesh zone name
                        mzName = obj.getDefaultMeshZoneName;
                        mzptr = emdlab_m3d_pmz(varargin{1},varargin{2});

                    else
                        error('Wrong input type.')
                    end

                case 4
                    % add a new mesh zone by direct passing name and connectivity list and points
                    mzName = obj.checkMeshZoneNonExistence(varargin{1});
                    if isnumeric(varargin{2}) && isnumeric(varargin{3})
                        mzptr = emdlab_m3d_pmz(varargin{2},varargin{3});
                    else
                        error('Wrong input type.')
                    end

                otherwise
                    error('Wrong inputs.');

            end

            % adding new mesh zone
            obj.mzs.(mzName) = mzptr;
            obj.mzs.(mzName).material = 'air';

            % changing states
            obj.isGlobalMeshGenerated = false;

        end

        function delete(obj)

            for mzName = obj.getMeshZoneNames
                delete(obj.mzs.(mzName));
            end

        end

        function y = get.isPRISM6(obj)
            y = strcmpi(obj.etype, 'PRISM6') || strcmpi(obj.etype, 'WEDGE6');
        end

        %% Solver Preparation Functions
        function ggmesh(obj, mzFlag)
            % generate global mesh

            if nargin<2, mzFlag = false; end

            % check states
            if obj.isGlobalMeshGenerated, return; end

            % generation of initial mesh
            Nn_tmp = 0;
            Ne_tmp = 0;
            mzNames = fieldnames(obj.mzs);

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                Nn_tmp = Nn_tmp + mzptr.Nn;
                Ne_tmp = Ne_tmp + mzptr.Ne;
            end

            % initialization of nodes and elements
            obj.nodes = zeros(Nn_tmp, 3);
            obj.cl = zeros(Ne_tmp, 6);
            obj.elements = zeros(Ne_tmp, 6);
            nindex = 0;
            eindex = 0;

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                % insertion of nodes
                obj.nodes(1 + nindex:mzptr.Nn + nindex, :) = mzptr.nodes;
                % insertion of elements
                obj.cl(1 + eindex:mzptr.Ne + eindex, :) = mzptr.cl + nindex;
                % specefying zone index
                obj.elements(1 + eindex:mzptr.Ne + eindex, 6) = i;
                mzptr.zi = i;
                nindex = nindex + mzptr.Nn;
                eindex = eindex + mzptr.Ne;
            end

            [obj.nodes, ~, ic] = uniquetol(obj.nodes, obj.gleps, 'ByRows', true);
            obj.cl = ic(obj.cl);
            % setting l2g
            nindex = 0;

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                mzptr.l2g = ic(nindex + 1:nindex + mzptr.Nn);
                nindex = nindex + mzptr.Nn;
            end

            obj.setData;
            obj.evalezi;

            if mzFlag

                for i = 1:obj.Nmzs
                    mzptr =  obj.mzs.(mzNames{i});
                    mzptr.props.cl = obj.cl(obj.ezi(:,mzptr.zi),:);
                end

            end

            % change states
            obj.isGlobalMeshGenerated = true;

            % settig element type
            obj.etype = 'PRISM6';

        end

        % evaluate element zone index
        function evalezi(obj)

            obj.ezi = false(obj.Ne, obj.Nmzs);

            for i = 1:obj.Nmzs
                obj.ezi(:, i) = obj.elements(:, 6) == i;
            end

        end

        % evaluate jacobian inverse transpose matrix
        function evalJIT(obj, mzFlag)

            if nargin<2, mzFlag = false; end

            % check states
            if obj.isJITEvaluated, return; end
            if obj.printFlag
                tic, disp('-------------------------------------------------------');
            end

            % prerequisite
            obj.ggmesh(mzFlag);

            % evaluation of jacobian inverse transpose using d1 element data
            % JIT is a [9 x Ne] matrix [J11;J21;J31;J12;J22;J32;J13;J23;J33]
            [obj.JIT, obj.gev] = emdlab_m3d_ttl4_evalJIT(obj.cl(:,1:4), obj.nodes);

            if mzFlag

                mzNames = fieldnames(obj.mzs);
                for i = 1:obj.Nmzs
                    mzptr = obj.mzs.(mzNames{i});
                    mzptr.props.JIT = obj.JIT(:,obj.ezi(:,mzptr.zi));
                    mzptr.props.gev = obj.gev(obj.ezi(:,mzptr.zi));
                end

            end

            % change states
            obj.isJITEvaluated = true;
            if obj.printFlag
                disp('Evaluation of JIT completed.');
                toc, disp('-------------------------------------------------------');
            end

        end

        %% Topological Functions
        function setData(obj)

            % getting the number of elements
            ne = obj.Ne;

            % prism facets
            f1 = [obj.cl(:, [1,2,3]), zeros(ne, 1)];
            f2 = [obj.cl(:, [4,6,5]), zeros(ne, 1)];
            f3 = obj.cl(:, [1,4,5,2]);
            f4 = obj.cl(:, [2,5,6,3]);
            f5 = obj.cl(:, [3,6,4,1]);

            % stack all faces
            allFacets = [f1; f2; f3; f4; f5];

            %             [~,idx] = min(allFacets,[],2);
            %             for i = 1:2*ne
            %                 allFacets(i,1:3) = circshift(allFacets(i,1:3), 1-idx(i));
            %             end
            %
            %             for i = (2*ne+1):5*ne
            %                 allFacets(i,:) = circshift(allFacets(i,:), 1-idx(i));
            %             end

            emdlab_mex_m3d_makeFacetsCanonical(allFacets);

            % canonical form for uniqueness
            sortedFacets = allFacets;
            s = false(size(allFacets,1),1);

            % triangles
            idx = (allFacets(:,2) > allFacets(:,3)) & (allFacets(:,4) == 0);
            sortedFacets(idx,2) = allFacets(idx,3);
            sortedFacets(idx,3) = allFacets(idx,2);
            s(idx) = true;

            % quadrilaterals
            idx = (allFacets(:,2) > allFacets(:,4)) & (allFacets(:,4) ~= 0);
            sortedFacets(idx,2) = allFacets(idx,4);
            sortedFacets(idx,4) = allFacets(idx,2);
            s(idx) = true;

            % unification of facets
            [obj.facets, ~, ic] = unique(sortedFacets, 'rows');
           
            % use negative sign for reverse facets
            ic(s) = -ic(s);

            % face index per element
            obj.elements(:,1:5) = reshape(ic, ne, 5);

            % calculate facet areas
            % initialize matrices
            obj.facetNormal = zeros(5*ne,3);
            obj.facetArea = zeros(5*ne,1);
            obj.fa = zeros(ne,5);
            obj.facetCenter = zeros(5*ne,3);

            % triangular facets
            idx = obj.facets(:,end) == 0;
            p12 = obj.nodes(obj.facets(idx,2),:) - obj.nodes(obj.facets(idx,1),:);
            p13 = obj.nodes(obj.facets(idx,3),:) - obj.nodes(obj.facets(idx,1),:);
            obj.facetNormal(idx,:) = cross(p12,p13);
            tmp = vecnorm(obj.facetNormal(idx,:),2,2);
            obj.facetNormal(idx,:) = obj.facetNormal(idx,:) ./ tmp;
            obj.facetArea(idx) = 0.5 * tmp;
            obj.fa(:,1:2) = obj.facetArea(abs(obj.elements(:,1:2)));
            obj.facetCenter(idx,:) = (obj.nodes(obj.facets(idx,1),:)+obj.nodes(obj.facets(idx,2),:)+...
                obj.nodes(obj.facets(idx,3),:))/3;

            % quadrilateral facets
            idx = ~idx;
            p12 = obj.nodes(obj.facets(idx,2),:) - obj.nodes(obj.facets(idx,1),:);
            p13 = obj.nodes(obj.facets(idx,3),:) - obj.nodes(obj.facets(idx,1),:);
            p14 = obj.nodes(obj.facets(idx,4),:) - obj.nodes(obj.facets(idx,1),:);
            obj.facetNormal(idx,:) = 0.5*(cross(p12,p13) + cross(p13,p14));
            obj.facetNormal(idx,:) = obj.facetNormal(idx,:) ./ vecnorm(obj.facetNormal(idx,:),2,2);
            obj.facetArea(idx) = 0.5*(vecnorm(cross(p12,p13),2,2) + vecnorm(cross(p13,p14),2,2));
            obj.fa(:,3:5) = obj.facetArea(abs(obj.elements(:,3:5)));
            obj.facetCenter(idx,:) = (obj.nodes(obj.facets(idx,1),:)+obj.nodes(obj.facets(idx,2),:)+...
                obj.nodes(obj.facets(idx,3),:)+obj.nodes(obj.facets(idx,4),:))/4;

            % facet center coordinates
            obj.xfc = obj.facetCenter(abs(obj.elements(:,1:5)),1);
            obj.yfc = obj.facetCenter(abs(obj.elements(:,1:5)),2);
            obj.zfc = obj.facetCenter(abs(obj.elements(:,1:5)),3);

            % calculate element centers
            obj.elementCenter = (obj.nodes(obj.cl(:,1),:) + obj.nodes(obj.cl(:,2),:) + obj.nodes(obj.cl(:,3),:) + obj.nodes(obj.cl(:,4),:) + ...
                obj.nodes(obj.cl(:,5),:) + obj.nodes(obj.cl(:,6),:))/6;

            % set nbs & facets
            obj.nbs = zeros(obj.Ne,5);
            obj.facets = [obj.facets, zeros(size(obj.facets,1), 6)];
            emdlab_m3d_pmdbc_evalfe(obj.facets, obj.elements, obj.nbs);

            % find boundary facets
            obj.bfacets = bitor(obj.facets(:,5) == 0, obj.facets(:,6) == 0);

            % last element of cl refers to mesh element type
            obj.cl(:,end+1) = 6;

            % last element of elements matrix refers to number of facets
            obj.elements(:,end+1) = 5;

            % last element of facets matrix refers to number of facet points
            obj.facets(idx,end+1) = 4;
            obj.facets(~idx,end) = 3;

            obj.NFN = 4;
            obj.NF = 5;
            obj.evalezi;
            obj.calculatePrismVolumes;

        end

        function calculatePrismVolumes(obj)
            %CALCULATEPRISMVOLUMES Compute volume of 6-node prism elements
            %
            % Inputs:
            %   obj.cl    : Ne x 6 connectivity
            %   obj.nodes : Nn x 3 nodal coordinates
            %
            % Prism node ordering assumed:
            %   bottom triangle: 1 2 3
            %   top triangle   : 4 5 6
            %
            % Shape functions:
            %   N1 = 0.5*(1-zeta)*(1-xi-eta)
            %   N2 = 0.5*(1-zeta)*xi
            %   N3 = 0.5*(1-zeta)*eta
            %   N4 = 0.5*(1+zeta)*(1-xi-eta)
            %   N5 = 0.5*(1+zeta)*xi
            %   N6 = 0.5*(1+zeta)*eta
            %
            % Reference domain:
            %   0 <= xi, eta, xi+eta <= 1
            %   -1 <= zeta <= 1

            Ne = size(obj.cl, 1);
            V  = zeros(Ne, 1);

            % 3-point quadrature on reference triangle
            % triangle area = 1/2, so weights = 1/6 each
            gpTri = [
                1/6, 1/6;
                2/3, 1/6;
                1/6, 2/3
                ];
            wTri = [1/6; 1/6; 1/6];

            % 2-point Gauss quadrature on [-1,1]
            g = 1 / sqrt(3);
            gpZ = [-g; g];
            wZ  = [1; 1];

            for e = 1:Ne
                Xe = obj.nodes(obj.cl(e,1:6), :);   % 6 x 3
                Ve = 0;

                for kz = 1:2
                    zeta = gpZ(kz);

                    for kt = 1:3
                        xi  = gpTri(kt,1);
                        eta = gpTri(kt,2);

                        % Derivatives of shape functions wrt xi
                        dN_dxi = 0.5 * [
                            -(1 - zeta)
                            (1 - zeta)
                            0
                            -(1 + zeta)
                            (1 + zeta)
                            0
                            ];

                        % Derivatives wrt eta
                        dN_deta = 0.5 * [
                            -(1 - zeta)
                            0
                            (1 - zeta)
                            -(1 + zeta)
                            0
                            (1 + zeta)
                            ];

                        % Derivatives wrt zeta
                        dN_dzeta = 0.5 * [
                            -(1 - xi - eta)
                            -xi
                            -eta
                            (1 - xi - eta)
                            xi
                            eta
                            ];

                        % Jacobian
                        J = [dN_dxi'; dN_deta'; dN_dzeta'] * Xe;   % 3x3

                        Ve = Ve + wTri(kt) * wZ(kz) * det(J);
                    end
                end

                % Positive volume even if element orientation is reversed
                V(e) = abs(Ve);
            end

            obj.gev = V(:)';
        end

        function idx = getFacetIndicesOnHalfPlane(obj, p0, nv)

            if ~obj.isdata
                obj.setData();
            end

            idx = [];

            for i = 1:size(obj.facets, 1)
                f = obj.facets(i, :);
                f = f(f > 0);
                xyz = obj.nodes(f, :);
                d = (xyz - p0) * nv(:);
                if all(abs(d) < obj.MINTOL)
                    idx(end + 1, 1) = i; %#ok<AGROW>
                end
            end
        end

        function idx = getFacetIndicesOnCylinder(obj, p0, p1, R, tol)
            %GETFACETINDICESONCYLINDER Returns indices of facets lying entirely on a cylinder's lateral surface.
            % A facet is on the cylinder if all of its nodes are on the cylinder.
            %
            % p0  - Point defining the start of the cylinder axis (1x3)
            % p1  - Point defining the end of the cylinder axis (1x3)
            % R   - Cylinder radius
            % tol - Tolerance (defaults to obj.gleps)

            % Handle default tolerance
            if nargin < 5 || isempty(tol)
                tol = obj.gleps;
            end

            % 1. Get indices of all nodes on the cylinder lateral surface
            node_idx = obj.getNodeIndicesOnCylinder(p0, p1, R, tol);

            % 2. Find facets whose all 4 nodes are in node_idx (using your ismember approach)
            mask = ismember(obj.facets(:, 1), node_idx) & ...
                ismember(obj.facets(:, 2), node_idx) & ...
                ismember(obj.facets(:, 3), node_idx) & ...
                ismember(obj.facets(:, 4), node_idx);

            % 3. Return facet indices
            idx = find(mask);
        end


        function varargout = showfb(obj)

            obj.ggmesh;
            [f,ax] = emdlab_r3d_geometry(1,0);

            cl_facets = obj.facets(obj.bfacets, 1:4);
            cl_facets(cl_facets == 0) = nan;
            patch('Faces', cl_facets, 'Vertices', obj.nodes, ...
                'FaceColor', 'b', 'EdgeColor', 'k', ...
                'FaceAlpha', 0.1, 'parent', ax);
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end
        
         %% Index Finding Functions
        function y = getfbf(obj)
            obj.ggmesh;
            y = find(obj.bfacets);
        end
        
        function showElement(obj, ie)

            if ~obj.ismesh
                obj.ggmesh();
            end

            if ie < 1 || ie > size(obj.elements, 1)
                error('Element index is out of range.');
            end

            cl = obj.elements(ie, :);

            F = [
                cl([1 2 3]) 0
                cl([4 5 6]) 0
                cl([1 2 5 4])
                cl([2 3 6 5])
                cl([3 1 4 6])
                ];

            [f, ax] = emdlab_r3d_geometry();
            f.Name = ['[Prism Element][', num2str(ie), ']'];
            hold(ax, 'on');

            FF = double(F);
            FF(FF == 0) = NaN;

            patch('Faces', FF, ...
                'Vertices', obj.nodes, ...
                'FaceColor', [0.8 0.9 1.0], ...
                'EdgeColor', 'k', ...
                'Parent', ax);

            axis(ax, 'equal');
            axis(ax, 'off');
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');
        end

        function varargout = showFacets(obj, varargin)

            [f,ax] = obj.showm;
            for i = 1:numel(varargin)
                tmp = obj.facets(varargin{i},[1,2,3,4]);
                tmp(tmp == 0) = nan;
                patch('Faces', tmp, 'Vertices', obj.nodes, ...
                    'FaceColor', 'r', 'EdgeColor', 'w', 'parent', ax, 'Marker', 'none', ...
                    'PickableParts', 'none', 'MarkerFaceColor', 'none', 'LineWidth', 1);
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function copyMirrorMeshZone(obj, idx, varargin)
            mz = copy(obj.mzs(idx));
            mz.mirror(varargin{:});
            mz.name = [mz.name, '_m'];
            obj.addMeshZone(mz, mz.name);
        end

        function cmmz(obj, idx, varargin)
            obj.copyMirrorMeshZone(idx, varargin{:});
        end

        function copyRotateMeshZone(obj, idx, varargin)
            mz = copy(obj.mzs(idx));
            mz.rotate(varargin{:});
            mz.name = [mz.name, '_r'];
            obj.addMeshZone(mz, mz.name);
        end

        function crmz(obj, idx, varargin)
            obj.copyRotateMeshZone(idx, varargin{:});
        end

        function copyShiftMeshZone(obj, idx, varargin)
            mz = copy(obj.mzs(idx));
            mz.shift(varargin{:});
            mz.name = [mz.name, '_s'];
            obj.addMeshZone(mz, mz.name);
        end

        function cshmz(obj, idx, varargin)
            obj.copyShiftMeshZone(idx, varargin{:});
        end

        function mirrorMeshZone(obj, idx, varargin)
            obj.mzs(idx).mirror(varargin{:});
            obj.ismesh = false;
            obj.isdata = false;
        end

        function mmz(obj, idx, varargin)
            obj.mirrorMeshZone(idx, varargin{:});
        end

        function rotateMeshZone(obj, idx, varargin)
            obj.mzs(idx).rotate(varargin{:});
            obj.ismesh = false;
            obj.isdata = false;
        end

        function rmz(obj, idx, varargin)
            obj.rotateMeshZone(idx, varargin{:});
        end

        function rotateMeshZones(obj, idxs, varargin)

            for i = 1:numel(idxs)
                obj.mzs(idxs(i)).rotate(varargin{:});
            end

            obj.ismesh = false;
            obj.isdata = false;
        end

        function shiftMeshZone(obj, idx, varargin)
            obj.mzs(idx).shift(varargin{:});
            obj.ismesh = false;
            obj.isdata = false;
        end

        function shmz(obj, idx, varargin)
            obj.shiftMeshZone(idx, varargin{:});
        end

        function shiftMeshZones(obj, idxs, varargin)

            for i = 1:numel(idxs)
                obj.mzs(idxs(i)).shift(varargin{:});
            end

            obj.ismesh = false;
            obj.isdata = false;
        end

        function joinMeshZones(obj, idxs, name)

            if nargin < 3
                name = "joined_mz";
            end

            idxs = idxs(:).';
            if isempty(idxs)
                return;
            end

            nodes = [];
            elements = [];
            offset = 0;

            for i = idxs
                mz = obj.mzs(i);
                nodes = [nodes; mz.nodes]; %#ok<AGROW>
                elements = [elements; mz.elements + offset]; %#ok<AGROW>
                offset = offset + mz.Nn;
            end

            mznew = emdlab_m3d_pmz(nodes, elements);
            mznew.name = char(name);

            obj.addMeshZone(mznew, name);
        end
    

    end

    methods (Access = private)

        function f = normalizeFace(obj, f)

            f = f(:).';
            nz = f(f > 0);

            if numel(nz) == 3
                nz = obj.minCyclic(nz);
                f = [nz, 0];
            elseif numel(nz) == 4
                nz = obj.minCyclic(nz);
                f = nz;
            else
                error('Facet must have 3 or 4 nodes.');
            end
        end

        function v = minCyclic(~, v)

            n = numel(v);
            cands = zeros(n, n);

            for i = 1:n
                cands(i, :) = circshift(v, [0, 1 - i]);
            end

            vr = fliplr(v);
            cands2 = zeros(n, n);
            for i = 1:n
                cands2(i, :) = circshift(vr, [0, 1 - i]);
            end

            allc = [cands; cands2];
            [~, idx] = min(allc * (max(allc(:)) + 1).^(n - 1:-1:0).');
            v = allc(idx, :);
        end
    end
end

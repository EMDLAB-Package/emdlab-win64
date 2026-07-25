% EMDLAB: Electrical Machines Design Laboratory
% common properties for all 3d mesh zone database classes

classdef emdlab_m3d_xmdb < handle & emdlab_mdb_cp

    properties

        % mesh nodes: [x,y]
        nodes (:,3) double;

        % mesh connectivity list: [n1, n2, n3, ...]
        cl (:,:) double;

        % mesh elements: [facet1, facet2, facet3, ...]
        elements (:,:) double;

        % Unique Mesh Facets
        facets (:,:) double;

        % Boundary Facets
        bfacets (:,1) logical;

        % edge length
        facetArea (:,1) double;
        fa (:,:) double;
        xfc (:,1) double;
        yfc (:,1) double;
        zfc (:,1) double;
        facetCenter (:,3) double;
        facetNormal (:,3) double;

        % center of elements
        elementCenter (:,3) double;

        % neighborhood elements
        nbs (:,:) double;

        % jacobian inverse transpose
        JIT (:,:) double;

        % global elements volume
        gev (1,:) double;

        % element zone index
        ezi (:,:) logical;

        % elements material index
        emi (:,:) logical;

        % auxiliary stored matricies
        mtcs (1,1) struct;

        % named selections
        facetNamedSelections (1,1) struct;

        % element type
        etype (1,:) char = '';

    end

    properties (Dependent = true)

        % Number of nodes
        Nn (1,1) double;

        % Number of elements
        Ne (1,1) double;

    end

    methods

        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end

        function y = get.Ne(obj)
            y = size(obj.cl, 1);
        end

        function ggmesh(~)
        end

        %% Visualization Functions
        function varargout = showm(obj, varargin)

            [f,ax] = emdlab_r3d_geometry(1,1);
            obj.ggmesh;
            mzNames = string(fieldnames(obj.mzs)');

            for mzName = mzNames
                mzptr = obj.mzs.(mzName);
                if isa(mzptr, 'emdlab_m3d_thmz')
                    plt = patch(ax,'Faces', obj.mzs.(mzName).facets(obj.mzs.(mzName).bfacets,1:3), ...
                        'Vertices', obj.mzs.(mzName).nodes, 'FaceColor', ...
                        obj.mzs.(mzName).color, 'EdgeColor', [0.2,0.2,0.2], ...
                        'FaceAlpha', 1, 'HitTest','on','PickableParts','visible');
                end
                plt.UserData.Tag = mzName;
                plt.UserData.c = obj.mzs.(mzName).color;
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showgg(obj, varargin)

            [f,ax] = emdlab_r3d_geometry(1,1);
            obj.ggmesh;
            mzNames = string(fieldnames(obj.mzs)');

            for mzName = mzNames
                mzptr = obj.mzs.(mzName);
                if isa(mzptr, 'emdlab_m3d_thmz')
                    plt = patch(ax,'Faces', obj.mzs.(mzName).facets(obj.mzs.(mzName).bfacets,1:3), ...
                        'Vertices', obj.mzs.(mzName).nodes, 'FaceColor', ...
                        obj.mzs.(mzName).color, 'EdgeColor', 'none', ...
                        'FaceAlpha', 1, 'HitTest','on','PickableParts','visible');
                end
                plt.UserData.Tag = mzName;
                plt.UserData.c = obj.mzs.(mzName).color;
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showg(obj, varargin)

            [f,ax] = emdlab_r3d_geometryNEW(1,0);
            obj.ggmesh;
            mzNames = string(fieldnames(obj.mzs)');

            for mzName = mzNames
                mzptr = obj.mzs.(mzName);
                if isa(mzptr, 'emdlab_m3d_thmz')
                    plt = patch(ax,'Faces', obj.mzs.(mzName).facets(obj.mzs.(mzName).bfacets,1:3), ...
                        'Vertices', obj.mzs.(mzName).nodes, 'FaceColor', ...
                        obj.mzs.(mzName).color, 'EdgeColor', 'none', ...
                        'FaceAlpha', 1, 'HitTest','on','PickableParts','visible');
                end
                plt.UserData.Tag = mzName;
                plt.UserData.c = obj.mzs.(mzName).color;
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showfb(obj)

            obj.ggmesh;
            [f,ax] = emdlab_r3d_geometry(1,0);

            patch('Faces', obj.facets(obj.bfacets, 1:3), 'Vertices', obj.nodes, ...
                'FaceColor', 'b', 'EdgeColor', 'k', ...
                'FaceAlpha', 0.1, 'parent', ax);
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = showwf(obj, varargin)
            % show wire frame mesh

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;

            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch('Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
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

        function varargout = showmzs(obj, varargin)
            % show mesh zones

            obj.showm;

        end

        function varargout = showmd(obj, varargin)
            % show mesh degree

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;

            mzNames = fieldnames(obj.mzs);
            ecolor = zeros(obj.Ne, 3);

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                ecolor(obj.ezi(:, i), 1) = mzptr.color(1);
                ecolor(obj.ezi(:, i), 2) = mzptr.color(2);
                ecolor(obj.ezi(:, i), 3) = mzptr.color(3);
            end

            % connectivity list index & point index
            switch obj.etype
                case {'T3', 'TL3'}
                    clIndex = 1:3;
                    pIndex = 1:3;
                case 'TL6'
                    clIndex = 1:3;
                    pIndex = 1:6;
                case 'TL10'
                    clIndex = 1:3;
                    pIndex = 1:10;
                case {'Q4', 'QL4'}
                    clIndex = 1:4;
                    pIndex = 1:4;
                otherwise
                    error('Wrong element type.');
            end

            patch('Faces', obj.cl(:, clIndex), 'Vertices', obj.nodes, ...
                'FaceColor', 'flat', 'FaceVertexCData', ecolor, ...
                'FaceAlpha', 0.5, 'EdgeColor', 'k', 'parent', ax);

            patch('Faces', obj.cl(:, pIndex), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'none', 'parent', ax, 'Marker', 'o', 'MarkerFaceColor', 'k');

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showNodes(obj, varargin)
            % show nodes on global mesh

            [f,ax] = obj.showm(varargin{1});
            if isnumeric(varargin{1})
                sIndex = 1;
            else
                sIndex = 2;
            end
            color = 'r';
            for i = sIndex:numel(varargin)
                patch('Faces', varargin{i}, 'Vertices', obj.nodes, ...
                    'EdgeColor', 'none', 'parent', ax, 'Marker', 'o', 'MarkerFaceColor', color,...
                    'HitTest', 'off','PickableParts','none');
                color = 'b';
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showEdges(obj, varargin)
            % show edges on global mesh

            [f,ax] = obj.showm;
            for i = 1:numel(varargin)
                patch('Faces', obj.edges(varargin{i},[1,2]), 'Vertices', obj.nodes, ...
                    'EdgeColor', 'r', 'parent', ax, 'Marker', 'o', 'MarkerFaceColor', 'r', ...
                    'LineWidth', 2,'HitTest', 'off','PickableParts','none');
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showElements(obj, varargin)
            % show elements on global mesh

            [f,ax] = obj.showm;

            % connectivity list index & point index
            switch obj.etype
                case {'T3', 'TL3'}
                    clIndex = 1:3;
                case 'TL6'
                    clIndex = 1:3;
                case 'TL10'
                    clIndex = 1:3;
                case {'Q4', 'QL4'}
                    clIndex = 1:4;
                otherwise
                    error('Wrong element type.');
            end

            for i = 1:numel(varargin)
                patch('Faces', obj.cl(varargin{i}, clIndex), 'Vertices', obj.nodes, ...
                    'FaceColor', 'r', 'EdgeColor', 'w', 'parent', ax, 'Marker', 'o', ...
                    'MarkerFaceColor', 'r', 'LineWidth', 2,...
                    'HitTest','off', 'PickableParts', 'none', 'handleVisibility', 'off');
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showFacets(obj, varargin)

            [f,ax] = obj.showm;
            for i = 1:numel(varargin)
                patch('Faces', obj.facets(varargin{i},[1,2,3]), 'Vertices', obj.nodes, ...
                    'FaceColor', 'r', 'EdgeColor', 'w', 'parent', ax, 'Marker', 'none', ...
                    'PickableParts', 'none', 'MarkerFaceColor', 'none', 'LineWidth', 1);
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showce(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;

            pts = [obj.getCenterOfElements; obj.getCenterOfEdges];
            cl1_tmp = repmat(1:obj.Ne,3,1);
            cl2_tmp = abs(obj.elements(:,1:3))' + obj.Ne;
            cl3_tmp = [cl1_tmp(:),cl2_tmp(:)];

            patch('Faces', obj.cl(:, [1,2,3]), 'Vertices', obj.nodes, ...
                'FaceColor', 'w', ...
                'FaceAlpha', 1, 'EdgeColor', 'k', 'linewidth', 1, 'parent', ax);

            patch('Faces', cl3_tmp, 'Vertices', pts, ...
                'FaceColor', 'none', 'EdgeColor', 'r', 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = shownes(obj)
            f = GraphicWindow;
            f.Name = '[Named edges of mesh]';
            h = guihandles(f);
            nes = fieldnames(obj.edgeNamedSelections);

            for i = 1:numel(nes)
                tmp = rand(1, 3);
                patch('Faces', obj.edges(obj.edgeNamedSelections.(nes{i}), 1:2), 'Vertices', ...
                    obj.nodes, 'FaceColor', ...
                    tmp, 'EdgeColor', tmp, 'linewidth', 1.5, ...
                    'FaceAlpha', 1, 'parent', h.va);
            end

            legend(h.va, nes);
            set(f, 'HandleVisibility', 'off', 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function showne(obj, neName)
            neName = obj.checkEdgeNamedSelectionExistence(neName);
            ah = setFigure('TMDBC: Named Edge');
            patch('Faces', obj.edges(obj.edgeNamedSelections.(neName), 1:2), 'Vertices', ...
                obj.nodes, 'FaceColor', ...
                'b', 'EdgeColor', 'b', 'linewidth', 1.5, ...
                'FaceAlpha', 1, 'parent', ah);
            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end

        function varargout = showvf(obj, Fx, Fy)
            f = GraphicWindow;
            f.Name = 'Field plot on center of elements';
            h = guihandles(f);
            h.va.NextPlot = 'add';
            c = obj.getCenterOfElements;
            color = zeros(obj.Ne, 3);
            mzNames = fieldnames(obj.mzs);

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                color(obj.ezi(:, mzptr.zi), :) = repmat(mzptr.color, mzptr.Ne, 1);
            end

            patch(h.va, 'faces', obj.cl(:, 1:3), 'vertices', obj.nodes, ...
                'facecolor', 'flat', 'FaceVertexCData', color, 'EdgeColor', 'w', 'facealpha', 0.5);
            quiver(h.va, c(:, 1), c(:, 2), Fx, Fy, 'color', 'k');
            axis(h.va, 'off', 'equal');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = showContact(obj, mz1Name, mz2Name)

            mz1Name = obj.checkMeshZoneExistence(mz1Name);
            mz2Name = obj.checkMeshZoneExistence(mz2Name);

            zi1 = obj.mzs.(mz1Name).zi;
            zi2 = obj.mzs.(mz2Name).zi;

            [f,ax] = emdlab_r3d_geometry(0,0);
            obj.ggmesh;

            idx = ((obj.facets(:, 4) == zi1) & (obj.facets(:, 5) == zi2)) | ...
                ((obj.facets(:, 4) == zi2) & (obj.facets(:, 5) == zi1));

            index = obj.facets(:, 4) ~= obj.facets(:, 5);
            patch('Faces', obj.facets(index, 1:3), 'Vertices', ...
                obj.nodes, 'FaceColor', ...
                'c', 'EdgeColor', 'none', ...
                'FaceAlpha', 0.2, 'parent', ax);

            patch('Faces', obj.facets(idx, 1:3), 'Vertices', ...
                obj.nodes, 'FaceColor', ...
                'y', 'EdgeColor', 'none', ...
                'FaceAlpha', 1, 'parent', ax);

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        %% Copy & Transform Functions
        function mirrorMeshZone(obj, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(mzName).mirror(varargin{:});
        end

        function mmz(varargin)
            mirrorMeshZone(varargin{:});
        end

        function copyMirrorMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getMirror(varargin{:});
        end

        function cmmz(varargin)
            copyMirrorMeshZone(varargin{:});
        end

        function rotateMeshZone(obj, mzName, varargin)
            mzName = char(mzName);
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(mzName).rotate(varargin{:});
        end

        function rotateMeshZones(obj, mzNames, varargin)
            for mzName = mzNames
                obj.rotateMeshZone(mzName, varargin{:});
            end
        end

        function rmz(varargin)
            rotateMeshZone(varargin{:});
        end

        function copyRotateMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getRotate(varargin{:});
        end

        function crmz(varargin)
            copyRotateMeshZone(varargin{:});
        end

        function shiftMeshZone(obj, mzName, varargin)
            mzName = erase(mzName, ' ');
            obj.checkMeshZoneExistence(mzName)
            obj.mzs.(mzName).shift(varargin{:});
        end

        function shiftMeshZones(obj, mzNames, shiftX, shiftY)
            for mzName = mzNames
                obj.shiftMeshZone(mzName, [shiftX, shiftY]);
            end
        end

        function shmz(varargin)
            shiftMeshZone(varargin{:});
        end

        function copyShiftMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getShift(varargin{:});
        end

        function cshmz(varargin)
            copyShiftMeshZone(varargin{:});
        end

        %% Boolean Functions
        function joinMeshZones(obj, nmzName, varargin)

            % find total number of mesh zones need to be joined
            xNmzs = 0;
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    xNmzs = xNmzs + 1;
                elseif isstring(varargin{i})
                    xNmzs = xNmzs + numel(varargin{i});
                else
                    error('Input type must be <char> or <string>.');
                end
            end

            if xNmzs < 2
                error('Minimum number mzs must be 2.');
            end

            mzNames = cell(1,xNmzs);
            index = 0;
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    index = index + 1;
                    mzNames{index} = obj.checkMeshZoneExistence(varargin{i});
                else
                    for j = 1:numel(varargin{i})
                        index = index + 1;
                        mzNames{index} = obj.checkMeshZoneExistence(char(varargin{i}(j)));
                    end
                end
            end

            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            Nn_tmp = zeros(1, xNmzs);
            Ne_tmp = zeros(1, xNmzs);

            for i = 1:xNmzs
                Nn_tmp(i) = obj.mzs.(mzNames{i}).Nn;
                Ne_tmp(i) = obj.mzs.(mzNames{i}).Ne;
            end

            n_nmz = zeros(sum(Nn_tmp), 2);
            e_nmz = zeros(sum(Ne_tmp), 4);
            n_tmp = 0;
            e_tmp = 0;

            for i = 1:xNmzs
                n_nmz(1 + n_tmp:n_tmp + Nn_tmp(i), :) = obj.mzs.(mzNames{i}).nodes;
                e_nmz(1 + e_tmp:e_tmp + Ne_tmp(i), :) = obj.mzs.(mzNames{i}).cl + n_tmp;
                n_tmp = n_tmp + Nn_tmp(i);
                e_tmp = e_tmp + Ne_tmp(i);
            end

            % jointing mzs
            [n_nmz, ~, ic] = uniquetol(n_nmz, obj.gleps, 'ByRows', true);
            e_nmz = ic(e_nmz);
            % adding new mz
            obj.mzs.(nmzName) = emdlab_m2d_qmz(e_nmz, n_nmz);
            obj.mzs.(nmzName).material = obj.mzs.(mzNames{1}).material;
            obj.mzs.(nmzName).color = obj.mzs.(mzNames{1}).color;
            % removing old mzs
            for i = 1:xNmzs
                obj.mzs = rmfield(obj.mzs, mzNames{i});
            end

        end

        function jmzs(varargin)
            joinMeshZones(varargin{:});
        end

        %% Index Finding Functions
        function y = getfbf(obj)
            % get indices of free boundary facets
            obj.ggmesh;
            y = find(obj.facets);
        end

        function y = getfbn(obj)
            y = obj.getfbf;
            y = obj.facets(y, 1:max(obj.elements(:,end)));
            y = unique(y(:));
            y = y(y~=0);
        end

        function idx = getNodeIndicesOnPlane(obj, p0, n, tol)
            %GETNODEINDICESONPLANE Indices of nodes lying on a plane.
            % Plane defined by point p0 and normal n.

            if nargin < 4 || isempty(tol)
                tol = obj.gleps;
            end

            % Ensure n is a column vector and p0 is a row vector for dot product
            n = n(:);
            p0 = p0(:).';

            % Compute dot product: Nodes * Normal - Point * Normal
            % This is the projection distance (scaled by ||n||)
            d = obj.nodes * n - (p0 * n);

            idx = find(abs(d) < tol);
        end

        function idx = getNodeIndicesOnHalfPlane(obj, p0, n, v_dir, tol)
            %GETNODEINDICESONHALFPLANE Indices of nodes on a half-plane.
            %
            % p0    - Point on the boundary line of the half-plane (1x3)
            % n     - Normal vector of the plane (1x3)
            % v_dir - Vector pointing into the active half of the plane (1x3)
            % tol   - Tolerance

            if nargin < 5 || isempty(tol)
                tol = obj.gleps;
            end

            % Force row vectors
            p0 = p0(:).';
            n = n(:).';
            v_dir = v_dir(:).';

            % Normalize vectors to make physical sense of tolerance
            n = n / norm(n);
            % Project v_dir to ensure it is orthogonal to the normal n
            v_dir = v_dir - (v_dir * n.') * n;
            v_dir = v_dir / norm(v_dir);

            % Vectors from p0 to all nodes (N x 3)
            V = obj.nodes - p0;

            % 1. Distance to the infinite plane (must be near 0)
            dist_to_plane = V * n.';

            % 2. Projection along the half-plane direction (must be >= 0)
            dist_along_half = V * v_dir.';

            % Nodes must lie on the plane AND on the positive side of the boundary line
            on_plane = abs(dist_to_plane) < tol;
            in_half = dist_along_half >= -tol; % allow tolerance at the boundary edge

            idx = find(on_plane & in_half);
        end

        function idx = getNodeIndicesOnCylinder(obj, p0, p1, h, tol)
            %GETNODEINDICESONCYLINDER Indices of nodes lying on a cylinder's lateral surface.
            % Cylinder axis goes from point p0 to point p1, with radius h.

            if nargin < 5 || isempty(tol)
                tol = obj.gleps;
            end

            % Ensure row vectors for vectorised calculations (1x3)
            p0 = p0(:).';
            p1 = p1(:).';

            % Axis vector and its length
            axis_vec = p1 - p0;
            L = norm(axis_vec);

            if L < 1e-12
                error('Points p0 and p1 must be distinct to define a cylinder axis.');
            end

            % Unit vector along the cylinder axis
            u = axis_vec / L;

            % Vectors from p0 to all nodes (N x 3)
            V = obj.nodes - p0;

            % Projection of vectors onto the axis (N x 1)
            d_axial = V * u.';

            % Squared distance to the axis: ||V||^2 - d_axial^2
            % sum(V.^2, 2) calculates the squared norm of each row
            d_perp_sq = sum(V.^2, 2) - d_axial.^2;

            % Ensure no negative values due to minor numerical precision limits
            d_perp = sqrt(max(d_perp_sq, 0));

            % Conditions:
            % 1. Radial distance must match radius h within tol
            % 2. Node must lie within the axial limits [0, L] of the cylinder
            on_lateral_surface = (abs(d_perp - h) < tol);
            within_bounds = (d_axial >= -tol) & (d_axial <= L + tol);

            idx = find(on_lateral_surface & within_bounds);
        end

        function idx = getFacetIndicesOnPlane(obj, p0, n, tol)
            %GETFACETINDICESONPLANE Returns indices of facets lying entirely on a plane.
            % A facet is on the plane if all of its nodes are on the plane.

            % Handle default tolerance
            if nargin < 4 || isempty(tol)
                tol = obj.gleps;
            end

            % Get indices of all nodes on the plane
            node_idx = obj.getNodeIndicesOnPlane(p0, n, tol);

            % Find facets whose all nodes are in node_idx

            mask = all(ismember(obj.facets(:,1:3), node_idx), 2);

            % Return facet indices
            idx = find(mask);
        end

        function idx = getFacetIndicesOnAnnulus(obj, p0, n, c0, rin, rout, tol)
            %GETFACETINDICESONANNULUS Returns indices of facets lying entirely on
            % an annular region of a plane.
            %
            % A facet is selected if all of its nodes:
            %   1) lie on the plane defined by point p0 and normal n
            %   2) lie between radii rin and rout from centre c0
            %
            % Inputs:
            %   p0   : 1x3 point on plane
            %   n    : 1x3 plane normal
            %   c0   : 1x3 centre of annulus (must lie on plane)
            %   rin  : inner radius
            %   rout : outer radius
            %   tol  : tolerance (optional)
            %
            % Output:
            %   idx  : indices of facets satisfying the condition

            % Default tolerance
            if nargin < 7 || isempty(tol)
                tol = obj.gleps;
            end

            % Normalise normal vector
            n = n(:).';
            n = n / norm(n);

            % Node coordinates
            % Assumes obj.nodes is Nx3
            X = obj.nodes;

            % Vector from plane point to nodes
            XP = X - p0;

            % Signed distance of each node from plane
            d = XP * n.';

            % Nodes on plane
            onPlane = abs(d) <= tol;

            % Vector from annulus centre to nodes
            XC = X - c0;

            % Remove normal component to get in-plane vector
            XC_plane = XC - (XC * n.') * n;

            % In-plane radial distance from annulus centre
            r = sqrt(sum(XC_plane.^2, 2));

            % Nodes inside annulus
            inAnnulus = (r >= rin - tol) & (r <= rout + tol);

            % Final node mask
            validNode = onPlane & inAnnulus;

            % Facet node indices (assuming triangular facets)
            F = obj.facets(:,1:3);

            % Keep facets whose all nodes are valid
            mask = all(validNode(F), 2);

            % Return facet indices
            idx = find(mask);
        end

        function idx = getEdgeIndicesOnContact(obj, mz1Name, mz2Name)

            mz1Name = obj.checkMeshZoneExistence(mz1Name);
            mz2Name = obj.checkMeshZoneExistence(mz2Name);

            mask = (ismember(obj.edges(:, 3), obj.mzs.(mz1Name).zi) & ...
                ismember(obj.edges(:, 4), obj.mzs.(mz2Name).zi)) | ...
                (ismember(obj.edges(:, 3), obj.mzs.(mz2Name).zi) & ...
                ismember(obj.edges(:, 4), obj.mzs.(mz1Name).zi));

            idx = find(mask);

        end

    end

end
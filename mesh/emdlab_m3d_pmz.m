% EMDLAB: Electrical Machines Design Laboratory
% Prism mesh zone (3D element)

classdef emdlab_m3d_pmz < handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_m3d_xmz

    methods
        %% Constructor and Destructor
        function obj = emdlab_m3d_pmz(cl, nodes)

            if nargin < 2, error('Not enough input arguments.'); end
            if nargin > 2, error('Too many input arguments.'); end

            obj.nodes = nodes;
            obj.cl = cl;
            obj.setData;

        end

        %% Solver Preparation Functions
        function evalVolumeOfElements(obj)

            % prism node ordering:
            % bottom triangle : 1 2 3
            % top triangle    : 4 5 6

            p1 = obj.nodes(obj.cl(:,1),:);
            p2 = obj.nodes(obj.cl(:,2),:);
            p3 = obj.nodes(obj.cl(:,3),:);
            p4 = obj.nodes(obj.cl(:,4),:);
            p5 = obj.nodes(obj.cl(:,5),:);
            p6 = obj.nodes(obj.cl(:,6),:);

            % decomposition of prism into three tetrahedra:
            % T1 = [1 2 3 4]
            % T2 = [2 3 5 4]
            % T3 = [3 5 6 4]

            v1 = evalTetraVolume(p1, p2, p3, p4);
            v2 = evalTetraVolume(p2, p3, p5, p4);
            v3 = evalTetraVolume(p3, p5, p6, p4);

            obj.ev = v1 + v2 + v3;

            function v = evalTetraVolume(p1, p2, p3, p4)
                v = abs(dot(p1-p4, cross(p2-p4, p3-p4, 2), 2)) / 6;
            end

        end

        function evalVolume(obj)

            obj.evalVolumeOfElements;
            obj.volume = sum(obj.ev);

        end

        %% Topological Functions
        function setdataForce(obj)
            obj.isDataSet = false;
            obj.setData;
        end

        function setData(obj)

            % check if data is already set
            if obj.isDataSet, return; end

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
            obj.elements = reshape(ic, ne, 5);
            
            % find boundary facets
            idx = emdlab_mex_findSignedPairs(obj.elements, size(obj.facets,1));
            obj.bfacets = idx ~= 3;

            % evaluation of volume
            obj.evalVolumeOfElements;
            obj.evalVolume;

            % update data set flag
            obj.isDataSet = true;

        end

        %% Mesh Visualization
        function varargout = showm(obj)

            [f,ax] = emdlab_r3d_geometry(1,0);
            f.Name = ['[Mesh Zone][', 'Nn = ', num2str(obj.Nn), '][Ne = ', num2str(obj.Ne), ']'];

            F = obj.facets(:,:);
            F = double(F);
            F(F == 0) = NaN;

            patch('Faces', F, 'Vertices', obj.nodes, 'FaceColor', ...
                obj.color, 'EdgeColor', 'k', 'parent', ax, 'facealpha', 1);

            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout == 2
                varargout{1} = f;
                varargout{2} = ax;
            elseif nargout > 2
                error('Too many output argument.');
            end

        end

        function showg(obj)
            obj.setData;

            f = GraphicWindow;
            h = guihandles(f);

            Ft = obj.getBoundaryTriFaces();

            patch('Faces', Ft, 'Vertices', obj.nodes, ...
                'FaceColor', obj.color, 'EdgeColor', 'none', 'parent', h.va);

            tr = triangulation(Ft, obj.nodes);

            try
                e = featureEdges(tr, pi/180);
            catch
                e = tr.edges;
            end

            patch('Faces', e, 'Vertices', tr.Points, ...
                'FaceColor', 'k', 'EdgeColor', 'k', 'parent', h.va);

            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end

        function varargout = showwf(obj)

            [f,ax] = emdlab_r3d_geometry(1,0);
            f.Name = ['[Mesh Zone][', 'Nn = ', num2str(obj.Nn), '][Ne = ', num2str(obj.Ne), ']'];

            F = obj.facets(obj.bfacets,:);
            F = double(F);
            F(F == 0) = NaN;

            patch('Faces', F, 'Vertices', obj.nodes, 'FaceColor', ...
                obj.color, 'EdgeColor', 'k', 'parent', ax, 'facealpha', 0.5);

            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout == 2
                varargout{1} = f;
                varargout{2} = ax;
            elseif nargout > 2
                error('Too many output argument.');
            end

        end

        %% Tools Functions
        function moveNodes(obj, MovTol)
            if nargin < 2
                MovTol = 1e-3;
            end

            obj.setData;

            % connectivity matrix for nodes
            Con = sparse(obj.edges(:,1), obj.edges(:,2), ...
                ones(size(obj.edges,1),1), obj.Nn, obj.Nn);
            Con = Con + Con';

            % loop for movements
            inodes = obj.getinodes;
            if isempty(inodes), return; end

            % weight matrix
            weight = spdiags(1./sum(Con(inodes,:),2), 0, numel(inodes), numel(inodes));

            for iter = 1:100
                % getting position of new nodes
                pnew = Con(inodes,:) * obj.nodes;
                pnew = weight * pnew;

                % evaluation of movements
                Mov = sqrt(sum((obj.nodes(inodes,:) - pnew).^2, 2));
                disp(sum(Mov));

                obj.nodes(inodes,:) = pnew;

                if all(Mov < MovTol)
                    disp(iter);
                    break
                end
            end

            obj.makeFalse_isDataSetted;
            obj.setData;
        end

        %% Index Finding
        function y = getfbf(obj)
            y = find(obj.bfacets);
        end

        function y = getfbn(obj)
            y = obj.getfbf;
            y = unique(y(:));
        end

        function y = getbnodes(obj)
            y = obj.facets(obj.bfacets,:);
            y = y(y > 0);
            y = unique(y(:));
        end

        function y = getinodes(obj)
            y = obj.getbnodes;
            y = setdiff((1:obj.Nn)', y);
        end

        function y = getNodeIndexOnPlane(obj, p0, n)
            y = find(abs(obj.nodes*n' - p0*n') < obj.gleps);
        end

        function y = getNodeIndexOnHalfPlane(obj, p0, p1, p2)
            y = obj.getNodeIndexOnPlane(p0, cross(p1, p2));
            tmp = obj.nodes(y,:) * p1' >= 0;
            y = y(tmp);
        end

        function y = getFacetIndexOnPlane(obj, varargin)
            yn = obj.getNodeIndexOnPlane(varargin{:});
            tf = false(size(obj.facets,1),1);

            for i = 1:size(obj.facets,1)
                face = obj.facets(i,:);
                face = face(face > 0);
                tf(i) = all(ismember(face, yn));
            end

            y = find(tf);
        end

        function y = getFacetIndexOnHalfPlane(obj, varargin)
            yn = obj.getNodeIndexOnHalfPlane(varargin{:});
            tf = false(size(obj.facets,1),1);

            for i = 1:size(obj.facets,1)
                face = obj.facets(i,:);
                face = face(face > 0);
                tf(i) = all(ismember(face, yn));
            end

            y = find(tf);
        end

        %% geometrical operations
        function y = getCenterOfElements(obj)
            y = (obj.nodes(obj.cl(:,1),:) + ...
                 obj.nodes(obj.cl(:,2),:) + ...
                 obj.nodes(obj.cl(:,3),:) + ...
                 obj.nodes(obj.cl(:,4),:) + ...
                 obj.nodes(obj.cl(:,5),:) + ...
                 obj.nodes(obj.cl(:,6),:)) / 6;
        end

        function y = getVolumeOfElements(obj)
            obj.setData;
            obj.evalVolumeOfElements;
            y = obj.ev;
        end

        function y = getVolume(obj)
            obj.setData;
            obj.evalVolumeOfElements;
            obj.evalVolume;
            y = obj.volume;
        end

        function y = getSurfaceArea(obj)

            bf = find(obj.bfacets);
            y = 0;

            for i = 1:numel(bf)
                face = obj.facets(bf(i),:);
                face = face(face > 0);

                if numel(face) == 3
                    p1 = obj.nodes(face(1),:);
                    p2 = obj.nodes(face(2),:);
                    p3 = obj.nodes(face(3),:);
                    y = y + 0.5 * norm(cross(p2 - p1, p3 - p1));
                else
                    p1 = obj.nodes(face(1),:);
                    p2 = obj.nodes(face(2),:);
                    p3 = obj.nodes(face(3),:);
                    p4 = obj.nodes(face(4),:);

                    a1 = 0.5 * norm(cross(p2 - p1, p3 - p1));
                    a2 = 0.5 * norm(cross(p3 - p1, p4 - p1));

                    y = y + a1 + a2;
                end
            end

        end

        %% tranforms and copy generations
        function mirror(obj, p0, p1)

            if nargin < 3
                p1 = p0;
                p0 = [0,0,0];
            end

            if ~isnumeric(p0) || ~isnumeric(p1)
                error('<p0> and <p1> must be numeric data.');
            elseif ~isequal(size(p0), [1,3]) || ~isequal(size(p1), [1,3])
                error('<p0> and <p1> must be 1x3 vectors.');
            end

            n = p1 - p0;
            n = n / norm(n);
            obj.nodes = obj.nodes - 2 * ((obj.nodes - p0) * n') * n;

        end

        function newObj = getMirror(obj, varargin)
            newObj = copy(obj);
            newObj.mirror(varargin{:});

            % reorder nodes to preserve prism orientation after mirroring
            newObj.cl = newObj.cl(:, [1,3,2,4,6,5]);

            newObj.makeFalse_isDataSetted;
            newObj.setData;
        end

        function newObj = getMirrorXY(obj)
            newObj = obj.getMirror([0,0,1]);
        end

        function newObj = getMirrorYZ(obj, varargin)
            newObj = obj.getMirror([1,0,0]);
        end

        function newObj = getMirrorZX(obj, varargin)
            newObj = obj.getMirror([0,1,0]);
        end

        function rotate(obj, rotAngle, p0, p1)

            if nargin < 4
                p1 = p0;
                p0 = [0,0,0];
            end

            if ~isnumeric(p0) || ~isnumeric(p1)
                error('<p0> and <p1> must be numeric data.');
            elseif ~isequal(size(p0), [1,3]) || ~isequal(size(p1), [1,3])
                error('<p0> and <p1> must be 1x3 vectors.');
            end

            u = p1 - p0;
            if norm(u) < eps
                error('Rotation axis is not defined. p0 and p1 must be distinct points.');
            end
            u = u / norm(u);

            V = obj.nodes - p0;

            cosT = cos(rotAngle);
            sinT = sin(rotAngle);

            dotUV   = V*u';
            crossUV = cross(repmat(u, size(V,1), 1), V, 2);

            V_rot = V*cosT + crossUV*sinT + dotUV*(1-cosT).*u;

            obj.nodes = V_rot + p0;
        end

        function newObj = getRotate(obj, varargin)
            newObj = copy(obj);
            newObj.rotate(varargin{:});
        end

        function newObj = getRotateX(obj, rotAngle)
            newObj = obj.getRotate(rotAngle, [1,0,0]);
        end

        function newObj = getRotateY(obj, rotAngle)
            newObj = obj.getRotate(rotAngle, [0,1,0]);
        end

        function newObj = getRotateZ(obj, rotAngle)
            newObj = obj.getRotate(rotAngle, [0,0,1]);
        end

        function shift(obj, p0, p1)

            if nargin < 3
                p1 = p0;
                p0 = [0,0,0];
            end

            if ~isnumeric(p0) || ~isnumeric(p1)
                error('<p0> and <p1> must be numeric data.');
            elseif ~isequal(size(p0), [1,3]) || ~isequal(size(p1), [1,3])
                error('<p0> and <p1> must be 1x3 vectors.');
            end

            d = p1 - p0;
            obj.nodes = obj.nodes + d;
        end

        function newObj = getShift(obj, varargin)
            newObj = copy(obj);
            newObj.shift(varargin{:});
        end

    end

    methods (Access = private)

        function makeFalse_isDataSetted(obj)
            obj.isDataSet = false;
            obj.isVolumeOfElementsEvaluated = false;
            obj.isMeshZoneVolumeEvaluated = false;
            obj.is_Wm_Evaluated = false;
            obj.is_Q_Evaluated = false;
        end

        function Ft = getBoundaryTriFaces(obj)

            bf = find(obj.bfacets);
            F = obj.facets(bf,:);

            triList = zeros(0,3);

            for i = 1:size(F,1)
                face = F(i,:);
                face = face(face > 0);

                if numel(face) == 3
                    triList(end+1,:) = face; %#ok<AGROW>
                else
                    triList(end+1,:) = face([1 2 3]); %#ok<AGROW>
                    triList(end+1,:) = face([1 3 4]); %#ok<AGROW>
                end
            end

            Ft = triList;
        end

    end

    methods (Static, Access = private)

        function [faceOut, sgn] = normalizeFace(faceIn)

            nodes = faceIn(faceIn > 0);
            n = numel(nodes);

            if ~(n == 3 || n == 4)
                error('Face must have 3 or 4 valid node indices.');
            end

            fwd = emdlab_m3d_pmz.minCyclic(nodes);
            rev = emdlab_m3d_pmz.minCyclic(fliplr(nodes));

            if emdlab_m3d_pmz.lexLessOrEqual(fwd, rev)
                base = fwd;
                sgn = 1;
            else
                base = rev;
                sgn = -1;
            end

            faceOut = zeros(1,4);
            faceOut(1:n) = base;
        end

        function vout = minCyclic(vin)

            n = numel(vin);
            cand = zeros(n, n);

            for k = 1:n
                cand(k,:) = circshift(vin, [0, 1-k]);
            end

            [~, idx] = sortrows(cand);
            vout = cand(idx(1),:);
        end

        function tf = lexLessOrEqual(a, b)

            tf = true;
            for i = 1:numel(a)
                if a(i) < b(i)
                    tf = true;
                    return;
                elseif a(i) > b(i)
                    tf = false;
                    return;
                end
            end
        end

        function edges = getUniqueEdgesFromFaces(facets)

            e = zeros(0,2);

            for i = 1:size(facets,1)
                face = facets(i,:);
                face = face(face > 0);

                if numel(face) == 3
                    tmp = [face([1 2]); face([2 3]); face([3 1])];
                elseif numel(face) == 4
                    tmp = [face([1 2]); face([2 3]); face([3 4]); face([4 1])];
                else
                    error('Face must have 3 or 4 nodes.');
                end

                tmp = sort(tmp, 2);
                e = [e; tmp]; %#ok<AGROW>
            end

            edges = unique(e, 'rows');
        end

    end

end

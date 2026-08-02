% EMDLAB: Electrical Machines Design Laboratory
% 2D triangular mesh zone

classdef emdlab_m2d_tmz < handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_m2d_xmz

    methods
        %% constructor and destructor
        function obj = emdlab_m2d_tmz(cl, nodes)

            if nargin < 2, error('Not enough input arguments.'); end
            if nargin > 2, error('Too many input arguments.'); end

            obj.nodes = nodes;
            obj.cl = cl;
            if ~isempty(cl)
                obj.setdata;
            end

        end

        %% FEM preparation
        function evalAreaOfElements(obj)
            %             if obj.is_ea_Evaluated, return; end
            v12 = obj.nodes(obj.cl(:, 2), :) - obj.nodes(obj.cl(:, 1), :);
            v13 = obj.nodes(obj.cl(:, 3), :) - obj.nodes(obj.cl(:, 1), :);
            obj.ea = 0.5 * (v12(:, 1) .* v13(:, 2) - v12(:, 2) .* v13(:, 1));
            % change states
            obj.isAreaOfElementsEvaluated = true;
        end

        function evalArea(obj)
            %             if obj.is_area_Evaluated, return; end
            obj.evalAreaOfElements;
            obj.area = sum(obj.ea);
            % change states
            obj.isMeshZoneAreaEvaluated = true;
        end

        %% topological functions
        % setting needed data
        function setdataForce(obj)
            obj.isDataSet = false;
            obj.setdata;
        end

        function setdata(obj)

            % check if already data is set
            if obj.isDataSet, return; end

            % getting the number of elements
            ne = obj.Ne;

            % triangle edges
            e1 = obj.cl(:,[1,2]);
            e2 = obj.cl(:,[2,3]);
            e3 = obj.cl(:,[3,1]);

            % stack all edges
            allEdges = [e1; e2; e3];
            sortedEdgets = allEdges;
            s = false(size(sortedEdgets,1),1);
            idx = allEdges(:,2) > allEdges(:,1);
            sortedEdgets(idx,1) = allEdges(idx,2);
            sortedEdgets(idx,2) = allEdges(idx,1);
            s(idx) = true;

            % unification of edges
            [obj.edges, ~, ic] = unique(sortedEdgets,'rows');

            % use negative sign for reverse edges
            ic(s) = -ic(s);

            % edge index per element
            obj.elements = reshape(ic, ne, 3);

            % find boundary edges
            idx = emdlab_mex_findSignedPairs(obj.elements, size(obj.edges,1));
            obj.bedges = idx ~= 3;

            % change states
            obj.isDataSet = true;

        end

        %% Named Selections
        function mzname = checkNodeNamedSelectionExistence(obj, mzname)
            mzname = rmspaces(mzname);

            if ~ isfield(obj.mzs, mzname)
                error('Specified node named selection does not exist.');
            end

        end

        function mzname = checkNodeNamedSelectionNonExistence(obj, mzname)
            mzname = rmspaces(mzname);

            if isfield(obj.mzs, mzname)
                error('Specified node named selection already exist.');
            end

        end

        function addNodeNamedSelection(obj, name, indices)
            name = obj.checkNodeNamedSelectionNonExistence(name);
            obj.nodeNamedSelections.(name) = indices;
        end

        %% Tools Functions
        function moveNodes(obj, iterMax, movTol)
            obj.setdata;

            if nargin < 2, iterMax = 100; end
            if nargin < 3, movTol = 1e-3; end

            % connectivity matrix for nodes
            Con = sparse(double(obj.edges(:, 1)), double(obj.edges(:, 2)), 1, obj.Nn, obj.Nn);
            Con = Con + Con';
            % loop for movments
            inodes = obj.getinodes;
            % weight matrix
            weight = diag(1 ./ sum(Con(inodes, :), 2));

            for iter = 1:iterMax
                % getting position of new nodes
                pnew = Con(inodes, :) * obj.nodes;
                pnew = weight * pnew;
                % evaluation of movments
                Mov = sqrt(sum((obj.nodes(inodes, :) - pnew).^2, 2));
                obj.nodes(inodes, :) = pnew;
                % check for movment tolerance
                if Mov < movTol
                    fprintf("Mesh smoothing #%d\n", iter);
                    break;
                end

            end

            % change states
            obj.isAreaOfElementsEvaluated = false;
            obj.isMeshZoneAreaEvaluated = false;
            obj.is_Q_Evaluated = false;
        end

        function moveNodes1(obj, iterMax, movTol)

            obj.setdata;
            if nargin < 2, iterMax = 100; end
            if nargin < 3, movTol = 1e-3; end

            % index of inner nodes
            idx = obj.getinodes;

            for iter = 1:iterMax

                % eij vectors
                eij = obj.nodes(obg.edges(:,2),:) - obj.nodes(obg.edges(:,1),:);

                % edge length
                el = vecnorm(eij,2,2);

                for i = idx
                end

                % getting position of new nodes
                pnew = Con(inodes, :) * obj.nodes;
                pnew = weight * pnew;
                % evaluation of movments
                Mov = sqrt(sum((obj.nodes(inodes, :) - pnew).^2, 2));
                obj.nodes(inodes, :) = pnew;
                % check for movment tolerance
                if Mov < movTol
                    fprintf("Mesh smoothing #%d\n", iter);
                    break;
                end

            end

            % change states
            obj.isAreaOfElementsEvaluated = false;
            obj.isMeshZoneAreaEvaluated = false;
            obj.is_Q_Evaluated = false;
        end

        function y = getbnodes(obj)
            % getting index of boundary nodes
            y = obj.edges(obj.bedges, :);
            y = unique(y(:));
        end

        function y = getinodes(obj)
            % getting index of inner nodes
            y = obj.getbnodes;
            y = setdiff((1:obj.Nn)', y);
        end

        function strefine(obj)

            % number of nodes in old mesh
            NnOld = obj.Nn;
            % nodes of new mesh
            obj.nodes = [obj.nodes; (obj.nodes(obj.edges(:, 1), :) + obj.nodes(obj.edges(:, 2), :)) / 2];
            % index of nodes on old edges
            index = abs(obj.elements);
            % new connctivity list
            obj.cl = [obj.cl(:, 1), index(:, [1, 3]) + NnOld
                obj.cl(:, 2), index(:, [2, 1]) + NnOld
                obj.cl(:, 3), index(:, [3, 2]) + NnOld
                index + NnOld];
            % setting data of new mesh
            obj.makeFalse_isDataSetted;
            obj.setdata;

        end

        %% tranforms and copy generations
        function mirror(obj, varargin)
            obj.nodes = ext_pmirror2(obj.nodes, varargin{:});
            obj.cl = obj.cl(:, [1, 3, 2]);
            obj.makeFalse_isDataSetted;
            obj.setdata;
        end

        function newObj = getMirror(obj, varargin)
            newObj = copy(obj);
            newObj.nodes = ext_pmirror2(newObj.nodes, varargin{:});
            newObj.cl = newObj.cl(:, [1, 3, 2]);
            newObj.makeFalse_isDataSetted;
            newObj.setdata;
        end

        function rotate(obj, varargin)
            if numel(varargin) == 1
                obj.nodes = ext_protate2(obj.nodes, varargin{:});
            else
                if length(varargin{2}) == 2
                    obj.nodes = ext_protate2(obj.nodes, varargin{:});
                else
                    obj.nodes = emdlab_g2d_rotatePoints(obj.nodes, varargin{:});
                end
            end
        end

        function newObj = getRotate(obj, varargin)
            newObj = copy(obj);
            newObj.nodes = ext_protate2(newObj.nodes, varargin{:});
        end

        function shift(obj, xShift, yShift)
            obj.nodes(:,1) = obj.nodes(:,1) + xShift;
            obj.nodes(:,2) = obj.nodes(:,2) + yShift;
        end

        function newObj = getShift(obj, varargin)
            newObj = copy(obj);
            newObj.nodes = ext_pshift2(newObj.nodes, varargin{:});
        end

        function evalQ(obj)
            %             if obj.is_Q_Evaluated, return; end
            obj.Q = sparse(double(obj.cl'), repmat(1:obj.Ne, 3, 1), ones(1, 3 * obj.Ne), ...
                obj.Nn, obj.Ne) * obj.ea;
            obj.Q = obj.Q' / 3 / obj.area;
            % change states
            obj.is_Q_Evaluated = true;
        end

        function evalWm(obj)
            %             if obj.is_Wm_Evaluated, return; end
            obj.Wm = sparse(double(obj.cl(:)), repmat((1:obj.Ne)', 3, 1), 1);
            % change states
            obj.is_Wm_Evaluated = true;
        end

        function smoothPlot(obj, value)
            [xdim, ydim] = size(value);

            if (xdim ~= obj.Ne) || (ydim ~= 1)
                error('Improper value, value must be a [Ne x 1] matrix.')
            end

            f = GraphicWindow(false);
            f.Name = 'Loss Density';
            h = guihandles(f);
            h.bg.Visible = 'off';
            h.ca.Visible = 'off';
            f.MenuBar = 'figure';
            f.Renderer = 'painters';
            f.Units = 'centimeters';
            f.Position = [10,10,13,13];

            patch('Faces', obj.cl, 'Vertices', obj.nodes, ...
                'FaceColor', 'interp', 'FaceVertexCdata', ...
                (obj.getWm * (value .* obj.getAreaOfElements))...
                ./ (obj.getWm * obj.getAreaOfElements), ...
                'EdgeColor', 'none', 'parent', h.va);

            patch('Faces', obj.edges(obj.bedges, :), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'parent', h.va);

            AddColorBar(f, min(value), max(value), 'loss density [W/Kg]', 8);
            set(f, 'Visible', 'on');

        end

    end

    methods (Access = private)

        function makeFalse_isDataSetted(obj)
            obj.isDataSet = false;
            obj.isAreaOfElementsEvaluated = false;
            obj.isMeshZoneAreaEvaluated = false;
            obj.is_Wm_Evaluated = false;
            obj.is_Q_Evaluated = false;
        end

    end

    %% Getters
    methods

        function y = getAreaOfElements(obj)
            obj.setdata;
            obj.evalAreaOfElements;
            y = obj.ea;
        end

        function y = getArea(obj)
            obj.setdata;
            obj.evalAreaOfElements;
            obj.evalArea;
            y = obj.area;
        end

        function y = getWm(obj)
            obj.setdata;
            obj.evalWm;
            y = obj.Wm;
        end

        function y = getQ(obj)
            obj.setdata;
            obj.evalAreaOfElements;
            obj.evalArea;
            obj.evalQ;
            y = obj.Q;
        end

        function y = getCenterOfElements(obj)
            % get center of elements
            y = (obj.nodes(obj.cl(:, 1), :) + ...
                obj.nodes(obj.cl(:, 2), :) + ...
                obj.nodes(obj.cl(:, 3), :)) / 3;
        end

        function y = getCenterOfEdges(obj)
            % get center of elements
            y = (obj.nodes(obj.edges(:, 1), :) + obj.nodes(obj.edges(:, 2), :)) / 2;
        end

        function y = getQuality(obj)
            % edges length
            el = sqrt(sum((obj.nodes(obj.edges(:, 1), :) - ...
                obj.nodes(obj.edges(:, 2), :)).^2, 2));
            b1 = el(abs(obj.elements(:, 1)));
            b2 = el(abs(obj.elements(:, 2)));
            b3 = el(abs(obj.elements(:, 3)));
            % mesh quality
            y = ((b1 + b2 - b3) .* (b1 - b2 + b3) .* (-b1 + b2 + b3)) ./ (b1 .* b2 .* b3);
            fprintf('Average Quality = %f\n', mean(y));
            fprintf('Minimum Quality = %f\n', min(y));
        end

        function y = getAspectRatio(obj)
            % edges length
            el = sqrt(sum((obj.nodes(obj.edges(:, 1), :) - ...
                obj.nodes(obj.edges(:, 2), :)).^2, 2));
            b1 = el(abs(obj.elements(:, 1)));
            b2 = el(abs(obj.elements(:, 2)));
            b3 = el(abs(obj.elements(:, 3)));
            y = max([b1, b2, b3], [], 2) ./ min([b1, b2, b3], [], 2);
        end

        function y = getEdgeLength(obj)
            % edges length
            y = sqrt(sum((obj.nodes(obj.edges(:, 1), :) - ...
                obj.nodes(obj.edges(:, 2), :)).^2, 2));
        end

        function y = getMaxEdgeLength(obj)
            % edges length
            el = sqrt(sum((obj.nodes(obj.edges(:, 1), :) - ...
                obj.nodes(obj.edges(:, 2), :)).^2, 2));
            b1 = el(abs(obj.elements(:, 1)));
            b2 = el(abs(obj.elements(:, 2)));
            b3 = el(abs(obj.elements(:, 3)));
            y = max([b1, b2, b3], [], 2);
        end

        function y = getCopy(obj)
            y = copy(obj);
        end

        function ttmz = getExtrude(obj, z, skewAngle)

            if iscolumn(z)
                z = z';
            end

            Nz = length(z);

            if nargin < 3
                z = repmat(z, obj.Nn, 1);
                ttmz = emdlab_m3d_thmz(tmzpc_getExtrude(obj.cl, obj.elements, ...
                    obj.Nn, Nz - 1), [repmat(obj.nodes, Nz, 1), z(:)]);
            else
                stepAngle = skewAngle * (pi / 180) / (Nz - 1);
                p = zeros(obj.Nn * Nz, 3);

                for i = 1:Nz
                    p((i - 1) * obj.Nn + 1:i * obj.Nn, 1:3) = ...
                        ext_protate3z([obj.nodes, repmat(z(i), obj.Nn, 1)], (i - 1) * stepAngle);
                end

                ttmz = TTMZPC(tmzpc_getExtrude(obj.cl, obj.elements, ...
                    obj.Nn, Nz - 1), p);
            end

        end

        function thmz = buildTetrahedralMeshByExtrusion(obj, z, skewAngle)

            if iscolumn(z)
                z = z';
            end

            Nz = length(z);
            if Nz < 2
                error('At least two z-levels are required to extrude.');
            end

            Ne2d = size(obj.cl, 1);
            Nlayers = Nz - 1;

            if nargin < 3 || isempty(skewAngle) || skewAngle == 0
                z3d = repmat(z, obj.Nn, 1);
                nodes3d = [repmat(obj.nodes, Nz, 1), z3d(:)];
            else
                stepAngle = skewAngle * (pi / 180) / (Nz - 1);
                nodes3d = zeros(obj.Nn * Nz, 3);

                for i = 1:Nz
                    p = [obj.nodes, repmat(z(i), obj.Nn, 1)];
                    ang = (i - 1) * stepAngle;

                    c = cos(ang);
                    s = sin(ang);

                    nodes3d((i - 1) * obj.Nn + 1:i * obj.Nn, :) = ...
                        [p(:,1) * c - p(:,2) * s, ...
                        p(:,1) * s + p(:,2) * c, ...
                        p(:,3)];
                end
            end

            % Each extruded triangle layer becomes one prism.
            % Each prism is split into 3 tetrahedra.
            cl3d = zeros(Ne2d * Nlayers * 3, 4);

            for k = 0:Nlayers-1
                e2d_idx = k * Ne2d + (1:Ne2d);
                nshift = k * obj.Nn;

                b1 = obj.cl(:,1) + nshift;
                b2 = obj.cl(:,2) + nshift;
                b3 = obj.cl(:,3) + nshift;

                t1 = b1 + obj.Nn;
                t2 = b2 + obj.Nn;
                t3 = b3 + obj.Nn;

                base = (k * Ne2d * 3);

                % Prism [b1 b2 b3 t1 t2 t3] split into 3 tetrahedra:
                % T1 = [b1 b2 b3 t1]
                % T2 = [b2 b3 t2 t1]
                % T3 = [b3 t2 t3 t1]
                cl3d(base + (1:Ne2d), :) = [b1, b2, b3, t1];
                cl3d(base + Ne2d + (1:Ne2d), :) = [b2, b3, t2, t1];
                cl3d(base + 2*Ne2d + (1:Ne2d), :) = [b3, t2, t3, t1];
            end

            thmz = emdlab_m3d_thmz(cl3d, nodes3d);

        end


        function mzptr = buildPrismMeshByExtrusion1(obj, z, skewAngle)

            if iscolumn(z)
                z = z';
            end

            Nz = length(z);
            if Nz < 2
                error('At least two z-levels are required to extrude.');
            end

            Ne2d = size(obj.cl, 1);
            Nlayers = Nz - 1;

            % Build prism connectivity
            cl3d = zeros(Ne2d * Nlayers, 6);
            for k = 0:Nlayers-1
                idx = k * Ne2d + (1:Ne2d);
                shift = k * obj.Nn;

                cl3d(idx,1:3) = obj.cl + shift;
                cl3d(idx,4:6) = obj.cl + shift + obj.Nn;
            end

            if nargin < 3 || isempty(skewAngle) || skewAngle == 0
                z3d = repmat(z, obj.Nn, 1);
                nodes3d = [repmat(obj.nodes, Nz, 1), z3d(:)];
            else
                stepAngle = skewAngle * (pi / 180) / (Nz - 1);
                nodes3d = zeros(obj.Nn * Nz, 3);

                for i = 1:Nz
                    p = [obj.nodes, repmat(z(i), obj.Nn, 1)];
                    ang = (i - 1) * stepAngle;

                    c = cos(ang);
                    s = sin(ang);

                    nodes3d((i-1)*obj.Nn + 1:i*obj.Nn, :) = ...
                        [p(:,1) * c - p(:,2) * s, ...
                        p(:,1) * s + p(:,2) * c, ...
                        p(:,3)];
                end
            end

            mzptr = emdlab_m3d_pmz(cl3d, nodes3d);
            mzptr.color = obj.color;
            mzptr.material = obj.material;

        end

 
function mzptr = buildPrismMeshByExtrusion(obj, z, skewAngle, draftAngle)
% buildPrismMeshByExtrusion
%
% Extrudes a 2D triangular mesh in the z-direction and optionally applies
% skew and draft during extrusion.
%
% Inputs:
%   z          : z-levels of the extrusion
%   skewAngle  : total skew angle in degrees
%                Positive value rotates the upper cross-section
%                counter-clockwise relative to the bottom one.
%
%   draftAngle : draft angle in degrees
%                Positive value expands the cross-section during extrusion.
%                Negative value contracts the cross-section.
%
% Output:
%   mzptr      : 3D prism mesh object
%
% Notes:
%   The draft is applied with respect to the centroid of the original
%   2D mesh. The radial displacement corresponding to the draft is
%
%       dr = dz * tan(draftAngle)
%
%   where dz is the distance from the bottom extrusion plane.
%
%   The skew and draft are applied simultaneously at every z-level.

    % -------------------------------------------------------------
    % Check and format z
    % -------------------------------------------------------------
    if iscolumn(z)
        z = z';
    end

    Nz = length(z);

    if Nz < 2
        error('At least two z-levels are required to extrude.');
    end

    % -------------------------------------------------------------
    % Default parameters
    % -------------------------------------------------------------
    if nargin < 3 || isempty(skewAngle)
        skewAngle = 0;
    end

    if nargin < 4 || isempty(draftAngle)
        draftAngle = 0;
    end

    % -------------------------------------------------------------
    % Basic mesh information
    % -------------------------------------------------------------
    Ne2d = size(obj.cl, 1);
    Nlayers = Nz - 1;

    % -------------------------------------------------------------
    % Build prism connectivity
    %
    % Each triangular element in the 2D mesh becomes a triangular
    % prism with 6 nodes:
    %
    %   bottom: 1 2 3
    %   top:    4 5 6
    % -------------------------------------------------------------
    cl3d = zeros(Ne2d * Nlayers, 6);

    for k = 0:Nlayers-1

        idx = k * Ne2d + (1:Ne2d);
        shift = k * obj.Nn;

        cl3d(idx,1:3) = obj.cl + shift;
        cl3d(idx,4:6) = obj.cl + shift + obj.Nn;

    end

    % -------------------------------------------------------------
    % 2D mesh centroid
    %
    % Draft is applied relative to this point.
    % -------------------------------------------------------------
    center = mean(obj.nodes, 1);

    % -------------------------------------------------------------
    % Convert angles to radians
    % -------------------------------------------------------------
    skewTotal = skewAngle * pi / 180;
    draftRad  = draftAngle * pi / 180;

    % -------------------------------------------------------------
    % Create 3D nodes
    % -------------------------------------------------------------
    nodes3d = zeros(obj.Nn * Nz, 3);

    % Total extrusion height
    H = z(end) - z(1);

    if H == 0
        error('The z-levels must span a non-zero extrusion height.');
    end

    for i = 1:Nz

        % Current z coordinate
        zi = z(i);

        % Distance from bottom plane
        dz = zi - z(1);

        % Normalized extrusion position
        t = dz / H;

        % ---------------------------------------------------------
        % Draft
        %
        % Radial displacement produced by the draft angle:
        %
        %     dr = dz * tan(draftAngle)
        %
        % To avoid a singularity at the centroid, the displacement
        % is applied by scaling the complete cross-section.
        %
        % The characteristic radius is taken as the maximum
        % distance of the original nodes from the centroid.
        % ---------------------------------------------------------
        if draftAngle == 0

            draftScale = 1;

        else

            r = sqrt( ...
                (obj.nodes(:,1) - center(1)).^2 + ...
                (obj.nodes(:,2) - center(2)).^2);

            R = max(r);

            if R == 0
                error('Cannot apply draft to a mesh with zero size.');
            end

            dr = dz * tan(draftRad);

            draftScale = 1 + dr / R;

            % Prevent inversion of the cross-section
            if draftScale <= 0
                error(['Draft angle is too large for the specified ', ...
                       'extrusion height. The cross-section would invert.']);
            end

        end

        % ---------------------------------------------------------
        % Apply draft
        % ---------------------------------------------------------
        x = center(1) + ...
            draftScale * (obj.nodes(:,1) - center(1));

        y = center(2) + ...
            draftScale * (obj.nodes(:,2) - center(2));

        % ---------------------------------------------------------
        % Apply skew
        %
        % The total skew angle is reached at the final z-level.
        % ---------------------------------------------------------
        ang = t * skewTotal;

        c = cos(ang);
        s = sin(ang);

        % Rotate around the mesh centroid rather than the global
        % origin. This is generally more appropriate for extrusion.
        xRel = x - center(1);
        yRel = y - center(2);

        xRot = center(1) + xRel * c - yRel * s;
        yRot = center(2) + xRel * s + yRel * c;

        % ---------------------------------------------------------
        % Store nodes
        % ---------------------------------------------------------
        idxNodes = (i-1)*obj.Nn + (1:obj.Nn);

        nodes3d(idxNodes,:) = [ ...
            xRot, ...
            yRot, ...
            repmat(zi, obj.Nn, 1)];

    end

    % -------------------------------------------------------------
    % Create prism mesh
    % -------------------------------------------------------------
    mzptr = emdlab_m3d_pmz(cl3d, nodes3d);

    % Copy properties
    mzptr.color = obj.color;
    mzptr.material = obj.material;

end



function [mzptr, lastLayerNodes] = buildPrismMeshByExtrusionZAxisDraft(obj, z, skewAngle, draftAngle)
% buildPrismMeshByExtrusionZAxisDraft
%
% Extrudes a 2D triangular mesh in the z-direction and optionally applies
% skew and radial draft with respect to the GLOBAL z-axis.
%
% Inputs:
%   z          : z-levels of the extrusion
%
%   skewAngle  : total skew angle in degrees
%                Positive value rotates the upper cross-section
%                counter-clockwise around the global z-axis.
%
%   draftAngle : draft angle in degrees
%
%                Positive value:
%                    cross-section SHRINKS toward the z-axis.
%
%                Negative value:
%                    cross-section EXPANDS away from the z-axis.
%
% Output:
%   mzptr          : 3D prism mesh object
%
%   lastLayerNodes : coordinates of the nodes at the final z-level
%                    after applying draft and skew.
%                    Size = [obj.Nn x 3]
%
% Description:
%
%   The draft is applied relative to the global z-axis:
%
%       x' = scale * x
%       y' = scale * y
%
%   The final radial displacement is defined by:
%
%       dr = H * tan(draftAngle)
%
%   where H is the total extrusion height.
%
%   Positive draftAngle shrinks the cross-section toward the z-axis.
%   Negative draftAngle expands the cross-section away from the z-axis.
%
%   Skew is applied as a rotation around the global z-axis.
%
%   Nodes located directly on the z-axis are handled naturally:
%
%       x = y = 0  -->  x' = y' = 0

    % -------------------------------------------------------------
    % Check and format z
    % -------------------------------------------------------------
    if iscolumn(z)
        z = z';
    end

    Nz = length(z);

    if Nz < 2
        error('At least two z-levels are required to extrude.');
    end

    % -------------------------------------------------------------
    % Default parameters
    % -------------------------------------------------------------
    if nargin < 3 || isempty(skewAngle)
        skewAngle = 0;
    end

    if nargin < 4 || isempty(draftAngle)
        draftAngle = 0;
    end

    % -------------------------------------------------------------
    % Basic mesh information
    % -------------------------------------------------------------
    Ne2d    = size(obj.cl, 1);
    Nlayers = Nz - 1;

    % -------------------------------------------------------------
    % Build prism connectivity
    %
    % Each triangular element becomes a 6-node prism:
    %
    %       bottom: 1 2 3
    %       top:    4 5 6
    % -------------------------------------------------------------
    cl3d = zeros(Ne2d * Nlayers, 6);

    for k = 0:Nlayers-1

        idx   = k * Ne2d + (1:Ne2d);
        shift = k * obj.Nn;

        cl3d(idx,1:3) = obj.cl + shift;
        cl3d(idx,4:6) = obj.cl + shift + obj.Nn;

    end

    % -------------------------------------------------------------
    % Convert angles to radians
    % -------------------------------------------------------------
    skewTotal = skewAngle * pi / 180;
    draftRad  = draftAngle * pi / 180;

    % -------------------------------------------------------------
    % Total extrusion height
    % -------------------------------------------------------------
    H = z(end) - z(1);

    if H == 0
        error('The z-levels must span a non-zero extrusion height.');
    end

    % -------------------------------------------------------------
    % Create 3D nodes
    % -------------------------------------------------------------
    nodes3d = zeros(obj.Nn * Nz, 3);

    % -------------------------------------------------------------
    % Calculate final radial scale
    %
    % R is the maximum radial distance from the global z-axis.
    % -------------------------------------------------------------
    if draftAngle ~= 0

        R = max(sqrt( ...
            obj.nodes(:,1).^2 + ...
            obj.nodes(:,2).^2));

        if R == 0
            error('The 2D mesh has zero radial size.');
        end

        % Radial displacement at the top layer
        dr = H * tan(draftRad);

        % New maximum radius
        Rnew = R - dr;

        % Prevent collapse or inversion
        if Rnew <= 0
            error(['Draft angle is too large for the specified ', ...
                   'extrusion height. The cross-section would ', ...
                   'collapse or invert.']);
        end

        % Scale at final layer
        scaleFinal = Rnew / R;

    else

        scaleFinal = 1;

    end

    % -------------------------------------------------------------
    % Generate each z-layer
    % -------------------------------------------------------------
    for i = 1:Nz

        % Current z-coordinate
        zi = z(i);

        % Distance from bottom plane
        dz = zi - z(1);

        % Normalized extrusion position
        t = dz / H;

        % ---------------------------------------------------------
        % Draft
        %
        % Linear interpolation from:
        %
        %       scale = 1
        %
        % at the bottom to:
        %
        %       scale = scaleFinal
        %
        % at the top.
        % ---------------------------------------------------------
        draftScale = 1 + t * (scaleFinal - 1);

        % ---------------------------------------------------------
        % Apply radial draft relative to global z-axis
        % ---------------------------------------------------------
        x = draftScale * obj.nodes(:,1);
        y = draftScale * obj.nodes(:,2);

        % ---------------------------------------------------------
        % Apply skew
        %
        % Total skewAngle is reached at the final layer.
        % ---------------------------------------------------------
        ang = t * skewTotal;

        c = cos(ang);
        s = sin(ang);

        xRot = c*x - s*y;
        yRot = s*x + c*y;

        % ---------------------------------------------------------
        % Store nodes
        % ---------------------------------------------------------
        idxNodes = (i-1)*obj.Nn + (1:obj.Nn);

        nodes3d(idxNodes,:) = [ ...
            xRot, ...
            yRot, ...
            repmat(zi, obj.Nn, 1)];

    end

    % -------------------------------------------------------------
    % Extract nodes of the LAST layer
    % -------------------------------------------------------------
    lastLayerNodes = nodes3d( ...
        (Nz-1)*obj.Nn + (1:obj.Nn), :);

    % -------------------------------------------------------------
    % Create prism mesh
    % -------------------------------------------------------------
    mzptr = emdlab_m3d_pmz(cl3d, nodes3d);

    % -------------------------------------------------------------
    % Copy properties
    % -------------------------------------------------------------
    mzptr.color    = obj.color;
    mzptr.material = obj.material;

end






        function ttmz = getRotateZ(obj, Angle, Nlayer)
            Angle = Angle * pi / 180;
            Angle = Angle / (Nlayer - 1);
            p = zeros(Nlayer * obj.Nn, 3);
            p(1:obj.Nn, 1:2) = obj.nodes;

            for i = 2:Nlayer
                p((1:obj.Nn) + (i - 1) * obj.Nn, :) = ...
                    ext_protate3z(p((1:obj.Nn) + (i - 2) * obj.Nn, :), Angle);
            end

            ttmz = TTMZPC(tmzpc_getExtrude(obj.cl, obj.elements, ...
                obj.Nn, Nlayer - 1), p);
        end

        function ttmz = getRotateY(obj, Angle, Nlayer)
            Angle = Angle * pi / 180;
            Angle = Angle / (Nlayer - 1);
            p = zeros(Nlayer * obj.Nn, 3);
            p(1:obj.Nn, 1:2) = obj.nodes;

            for i = 2:Nlayer
                p((1:obj.Nn) + (i - 1) * obj.Nn, :) = ...
                    ext_protate3y(p((1:obj.Nn) + (i - 2) * obj.Nn, :), Angle);
            end

            ttmz = TTMZPC(tmzpc_getExtrude(obj.cl, obj.elements, ...
                obj.Nn, Nlayer - 1), p);
        end

        function ttmz = getRotateX(obj, Angle, Nlayer)
            Angle = Angle * pi / 180;
            Angle = Angle / (Nlayer - 1);
            p = zeros(Nlayer * obj.Nn, 3);
            p(1:obj.Nn, 1:2) = obj.nodes;

            for i = 2:Nlayer
                p((1:obj.Nn) + (i - 1) * obj.Nn, :) = ...
                    ext_protate3x(p((1:obj.Nn) + (i - 2) * obj.Nn, :), Angle);
            end

            ttmz = TTMZPC(tmzpc_getExtrude(obj.cl, obj.elements, ...
                obj.Nn, Nlayer - 1), p);
        end

    end

end

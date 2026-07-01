% EMDLAB: Electrical Machines Design Laboratory
% 2D triangular mesh zone

classdef emdlab_m2d_tmz < handle & emdlab_g2d_constants & matlab.mixin.Copyable
    
    properties
        
        % mesh nodes
        nodes (:,2) double;

        % mesh connectivity list
        cl (:,3) double;

        % mesh elements: [edge1, edge2, edge3]
        elements (:,3) double;

        % unique edges: [node1, node2]
        edges (:,2) double;

        % list of boundary edges
        bedges (:,1);

        % Named Selections
        nodeNamedSelections (1,1) struct;
        edgeNamedSelections (1,1) struct;
        
    end
    
    properties
        
        % zone index
        zi (1,1) double;

        % local to global node index
        l2g (:,1) double;

        % material of zone
        material char = 'air';

        % mesh zone color
        color = 'c';

        % surface color transparency
        transparency (1,1) double = 1;
        
        % mesh zone properties: differs in differents solvers
        props (1,1) struct;

        % flags
        isMoving = false;
        
    end
    
    properties (Access = private)
        
        % a vector containing area of elements
        ea (:,1) double;

        % mesh zone area
        area (1,1) double;

        % Q
        Q (1,:) double;

        % Weight matrix
        Wm

        % states
        isDataSetted (1,1) logical = false;
        is_ea_Evaluated (1,1) logical = false;
        is_area_Evaluated (1,1) logical = false;
        is_Wm_Evaluated (1,1) logical = false;
        is_Q_Evaluated (1,1) logical = false;
        
    end
    
    properties (Dependent = true)
        
        % number of mesh zone nodes
        Nn (1,1) double;

        % number of mesh zone elements
        Ne (1,1) double;
        
    end
    
    methods        
        %% constructor and destructor
        function obj = emdlab_m2d_tmz(cl, nodes)

            if nargin < 2, error('Not enough input arguments.'); end
            if nargin > 2, error('Too many input arguments.'); end

            obj.nodes = nodes;
            obj.cl = cl;
            obj.setdata;

        end
        
        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end
        
        function y = get.Ne(obj)
            y = size(obj.cl, 1);
        end
              
        %% FEM preparation
        function evalAreaOfElements(obj)
%             if obj.is_ea_Evaluated, return; end
            v12 = obj.nodes(obj.cl(:, 2), :) - obj.nodes(obj.cl(:, 1), :);
            v13 = obj.nodes(obj.cl(:, 3), :) - obj.nodes(obj.cl(:, 1), :);
            obj.ea = 0.5 * (v12(:, 1) .* v13(:, 2) - v12(:, 2) .* v13(:, 1));
            % change states
            obj.is_ea_Evaluated = true;
        end
        
        function evalArea(obj)
%             if obj.is_area_Evaluated, return; end
            obj.evalAreaOfElements;
            obj.area = sum(obj.ea);
            % change states
            obj.is_area_Evaluated = true;
        end
        
        %% topological functions
        % setting needed data
        function setdataForce(obj)
            obj.isDataSetted = false;
            obj.setdata;
        end

        function setdata(obj)

            % check if already data is set
            if obj.isDataSetted, return; end

            % first edge of each triangle
            e1 = obj.cl(:, [1, 2]);
            % second edge of each triangle
            e2 = obj.cl(:, [2, 3]);
            % third edge of each triangle
            e3 = obj.cl(:, [3, 1]);

            % sorting for lower index
            [e1, s1] = sort(e1, 2);
            [e2, s2] = sort(e2, 2);
            [e3, s3] = sort(e3, 2);

            % specefying changed edge index
            s1 = s1(:, 1) == 2;
            s2 = s2(:, 1) == 2;
            s3 = s3(:, 1) == 2;
            
            % unification of edges
            [obj.edges, ~, ic] = unique([e1; e2; e3], 'rows');

            % getting number of elements
            ne = obj.Ne;

            % getting index of edge corresponding to each elements
            e1 = ic(1:ne);
            e2 = ic(1 + ne:2 * ne);
            e3 = ic(1 + 2 * ne:3 * ne);

            % specefying boundary edges
            obj.bedges = sparse([e1, e2, e3], ones(3 * ne, 1), ones(3 * ne, 1));
            obj.bedges = full(obj.bedges == 1);

            % specefying trace direction
            e1(s1) = -e1(s1);
            e2(s2) = -e2(s2);
            e3(s3) = -e3(s3);

            % element matrix
            obj.elements = [e1, e2, e3];

            % change states
            obj.isDataSetted = true;

        end
        
        %% visiualization functions
        function varargout = showm(obj)

            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = ['[Mesh Zone][', 'Nn = ', num2str(obj.Nn), '][Ne = ', num2str(obj.Ne), ']'];
            
            patch('Faces', obj.cl(:, 1:3), 'Vertices', obj.nodes, 'FaceColor', ...
                obj.color, 'EdgeColor', 'k', 'parent', ax);
            
            zoom on;
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
        
        function f = showwf(obj)
            f = GraphicWindow();
            f.Name = ['[Mesh Zone Boundary Edges][', 'Nbe = ', num2str(length(obj.bedges)), ']'];
            h = guihandles(f);
            patch('Faces', obj.edges(obj.bedges, [1,2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'parent', h.va);
            set(f, 'HandleVisibility', 'off', 'Visible', 'on');
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
        function moveNodes(obj, MovTol)
            obj.setdata;
            
            if nargin < 2
                MovTol = 1e-3;
            end
            
            % connectivity matrix for nodes
            Con = sparse(double(obj.edges(:, 1)), double(obj.edges(:, 2)), 1, obj.Nn, obj.Nn);
            Con = Con + Con';
            % loop for movments
            inodes = obj.getinodes;
            % weight matrix
            weight = diag(1 ./ sum(Con(inodes, :), 2));
            
            for iter = 1:100
                % getting position of new nodes
                pnew = Con(inodes, :) * obj.nodes;
                pnew = weight * pnew;
                % evaluation of movments
                Mov = sqrt(sum((obj.nodes(inodes, :) - pnew).^2, 2));
                obj.nodes(inodes, :) = pnew;
                % check for movment tolerance
                if Mov < MovTol
                    disp(iter);
                    break;
                end
                
            end
            
            % change states
            obj.is_ea_Evaluated = false;
            obj.is_area_Evaluated = false;
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
            obj.isDataSetted = false;
            obj.is_ea_Evaluated = false;
            obj.is_area_Evaluated = false;
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
                ttmz = TTMZPC(tmzpc_getExtrude(obj.cl, obj.elements, ...
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

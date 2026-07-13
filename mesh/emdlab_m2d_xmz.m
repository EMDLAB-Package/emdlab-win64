% EMDLAB: Electrical Machines Design Laboratory
% common properties for all mesh 2d mesh zone classes

classdef emdlab_m2d_xmz < handle

    properties

        % mesh nodes
        nodes (:,2) double;

        % mesh connectivity list
        % [n1, n2, n3] -> for triangular mesh
        % [n1, n2, n3, n4] -> for quadrilateral mesh
        % [n1, n2, n3, 0 or n4, nEdges] -> for qt mesh
        % [n1, n2, n3, ..., nEdges] -> for arbitrary mesh
        cl (:,:) double;

        % mesh elements
        % mesh elements: [edge1, edge2, edge3] -> for triangular mesh
        % mesh elements: [edge1, edge2, edge3, edge4] -> for quadrilateral mesh
        % mesh elements: [edge1, edge2, edge3, 0 or edge4, nEdges] -> for qt mesh
        % mesh elements: [edge1, edge2, edge3, ..., nEdges] -> for arbitrary mesh
        elements (:,:) double;

        % unique edges: [node1, node2]
        edges (:,2) double;

        % list of boundary edges
        bedges (:,1);

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

        % mesh zone orientation
        orientation = 'global';

        % a vector containing area of elements
        ea (:,1) double;

        % mesh zone area
        area (1,1) double;

        % states
        isDataSet (1,1) logical = false;
        isAreaOfElementsEvaluated (1,1) logical = false;
        isMeshZoneAreaEvaluated (1,1) logical = false;

        % Q
        Q (1,:) double;

        % Weight matrix
        Wm

        % states
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
        %% Get Dependents
        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end

        function y = get.Ne(obj)
            y = size(obj.cl, 1);
        end

        %% Flag manipulations
        function clearSetDataFlag(obj)

            obj.isDataSet = false;
            obj.isAreaOfElementsEvaluated = false;
            obj.isMeshZoneAreaEvaluated = false;
            obj.is_Wm_Evaluated = false;
            obj.is_Q_Evaluated = false;

        end

        %% Visiualization Functions
        function varargout = showm(obj)

            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = "Nn = " + string(obj.Nn) + ", Ne" + string(obj.Ne);

            if ismember(size(obj.cl,2), [3,4])
                patch('Faces', obj.cl, 'Vertices', obj.nodes, 'FaceColor', ...
                    obj.color, 'EdgeColor', 'k', 'parent', ax);
            end

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

        function varargout = showwf(obj)

            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = "Nn = " + string(obj.Nn) + ", Ne" + string(obj.Ne);

            patch('Faces', obj.edges(obj.bedges,1:2), 'Vertices', obj.nodes, 'FaceColor', ...
                obj.color, 'EdgeColor', 'k', 'parent', ax, 'linewidth', 1);

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

        %% Tranform & Copy
        function mirror(obj, varargin)
            obj.nodes = ext_pmirror2(obj.nodes, varargin{:});
            if size(obj.cl,2) == 3
                obj.cl = obj.cl(:, [1, 3, 2]);
            elseif size(obj.cl,2) == 4
                obj.cl = obj.cl(:, [1, 4, 3, 2]);
            else
            end
            obj.clearSetDataFlag;
            obj.setData;
        end

        function newObj = getMirror(obj, varargin)
            newObj = copy(obj);
            newObj.nodes = ext_pmirror2(newObj.nodes, varargin{:});
            if size(obj.cl,2) == 3
                obj.cl = obj.cl(:, [1, 3, 2]);
            elseif size(obj.cl,2) == 4
                obj.cl = obj.cl(:, [1, 4, 3, 2]);
            else
            end
            newObj.clearSetDataFlag;
            newObj.setData;
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

    end

end
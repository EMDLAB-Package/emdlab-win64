% EMDLAB: Electrical Machines Design Laboratory
% Quadrilateral-Triangular mesh zone (2D element)

classdef emdlab_m2d_qtmz <  handle & emdlab_g2d_constants & matlab.mixin.Copyable

    properties (SetAccess = private)

        % mesh nodes
        nodes (:,2) double;

        % quadrilateral mesh connectivity list
        qcl (:,4) double;

        % triangular mesh connectivity list
        tcl (:,3) double;

        % quadrilateral mesh elements
        qelements (:,4) double;

        % triangular mesh elements
        telements (:,3) double;

        % unique edges
        edges (:,2) double;

        % list of boundary edges
        bedges (:,1) double;

    end

    properties

        % zone index
        zi (1,1) double;

        % local to global node index
        l2g (:,1) double;

        % material of zone
        material char = 'air';

        % surface color transparency
        transparency (1,1) double = 1;

        % mesh zone color
        color = 'c';

        % mesh zone properties: differs in differents solvers
        props (1,1) struct;

        % flags
        isMoving = false;

        orientation = 'global';

    end

    properties (Access = private)

        % a vector containing area of elements
        ea (:,1) double;

        % mesh zone area
        area (1,1) double;

        % states
        isDataSet (1,1) logical = false;
        isAreaOfElementsEvaluated (1,1) logical = false;
        isMeshZoneAreaEvaluated (1,1) logical = false;

    end

    properties (Dependent = true)

        % number of mesh zone nodes
        Nn (1,1) double;

        % number of quadrilateral elements
        Nq (1,1) double;

        % number of triangular elements
        Nt (1,1) double;

    end

    methods
        %% Constructor and Destructor
        function obj = emdlab_m2d_qtmz(qcl, tcl, nodes)

            if nargin < 3, error('Not enough input arguments.'); end
            if nargin > 3, error('Too many input arguments.'); end

            obj.nodes = nodes;
            obj.qcl = qcl;
            obj.tcl = tcl;
            obj = obj.setData;

        end

        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end

        function y = get.Nq(obj)
            y = size(obj.qcl, 1);
        end

        function y = get.Nt(obj)
            y = size(obj.tcl, 1);
        end

        %% Topological Functions
        % setting needed data
        function obj = setData(obj)

            % check if already data is set
            if obj.isDataSet, return; end

            % first edge of each quadrilateral
            e1 = obj.qcl(:,[1,2]);
            % second edge of each quadrilateral
            e2 = obj.qcl(:,[2,3]);
            % third edge of each quadrilateral
            e3 = obj.qcl(:,[3,4]);
            % forth edge of each quadrilateral
            e4 = obj.qcl(:,[4,1]);

            % sorting for lower index
            [e1,s1] = sort(e1, 2);
            [e2,s2] = sort(e2, 2);
            [e3,s3] = sort(e3, 2);
            [e4,s4] = sort(e4, 2);

            % specefying changed edge index
            s1 = s1(:,1) == 2;
            s2 = s2(:,1) == 2;
            s3 = s3(:,1) == 2;
            s4 = s4(:,1) == 2;

            % unification of edges
            [obj.edges, ~, ic] = unique([e1; e2; e3; e4],'rows');

            % getting number of elements
            ne = obj.Nq;

            % getting index of edge corresponding to each elements
            e1 = ic(1:ne);
            e2 = ic(1+ne:2*ne);
            e3 = ic(1+2*ne:3*ne);
            e4 = ic(1+3*ne:4*ne);

            % specefying boundary edges
            obj.bedges = sparse([e1,e2,e3,e4],ones(4*ne,1),ones(4*ne,1));
            obj.bedges = full(obj.bedges==1);

            % specefying trace direction
            e1(s1) = -e1(s1);
            e2(s2) = -e2(s2);
            e3(s3) = -e3(s3);
            e4(s4) = -e4(s4);

            % element matrix
            obj.qelements = [e1,e2,e3,e4];

            % change states
            obj.isDataSet = true;

        end

        %% Mesh Visiualization
        function varargout = showm(obj)

            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = "Nn = " + string(obj.Nn) + ", Nq = " + string(obj.Nq) + ", Nt = " + string(obj.Nt);

            patch('Faces',obj.qcl,'Vertices',obj.nodes,'FaceColor',...
                obj.color,'EdgeColor','k','parent',ax);

            patch('Faces',obj.tcl,'Vertices',obj.nodes,'FaceColor',...
                obj.color,'EdgeColor','k','parent',ax);

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

            patch('Faces',obj.edges(find(obj.bedges),:),'Vertices',obj.nodes,...
                'FaceColor','none','EdgeColor','k', 'parent', ax);
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

        %% Tools Functions
        function obj = moveNodes(obj,MovTol)
            if nargin<2
                MovTol = 1e-3;
            end
            % connectivity matrix for nodes
            Con = sparse(obj.edges(:,1),obj.edges(:,2),...
                ones(size(obj.edges,1),1),obj.Nn,obj.Nn);
            Con = Con + Con';
            % loop for movments
            inodes = obj.getInnerNodes;
            % weight matrix
            weight = diag(1./sum(Con(inodes,:),2));
            for iter = 1:100
                % getting position of new nodes
                pnew = Con(inodes,:)*obj.nodes;
                pnew = weight*pnew;
                % evaluation of movments
                Mov = sqrt(sum((obj.nodes(inodes,:)-pnew).^2,2));
                disp(sum(Mov));
                obj.nodes(inodes,:) = pnew;
                % check for movment tolerance
                if Mov < MovTol
                    disp(iter);
                    break
                end
            end
        end

        function y = getBundaryNodes(obj)

            % getting index of boundary nodes
            y = obj.edges(obj.bedges,:);
            y = unique(y(:));

        end

        function y = getInnerNodes(obj)

            % getting index of inner nodes
            y = obj.getBundaryNodes;
            y = setdiff((1:obj.Nn)',y);

        end

        function strefine(obj)

            % number of nodes and elements in old mesh
            NnOld = obj.Nn;
            NeOld = obj.Nq;
            NedOld = size(obj.edges, 1);
            % nodes of new mesh
            obj.nodes = [obj.nodes;...
                (obj.nodes(obj.edges(:,1),:)+obj.nodes(obj.edges(:,2),:))/2;
                (obj.nodes(obj.qcl(:,1),:)+obj.nodes(obj.qcl(:,2),:)+obj.nodes(obj.qcl(:,3),:)+obj.nodes(obj.qcl(:,4),:))/4];
            % index of nodes on old edges
            index = [abs(obj.qelements),(1:NeOld)'];
            % new connctivity list
            obj.qcl = [obj.qcl(:,1),index(:,1)+NnOld,index(:,5)+NnOld+NedOld,index(:,4)+NnOld
                obj.qcl(:,2),index(:,2)+NnOld,index(:,5)+NnOld+NedOld,index(:,1)+NnOld
                obj.qcl(:,3),index(:,3)+NnOld,index(:,5)+NnOld+NedOld,index(:,2)+NnOld
                obj.qcl(:,4),index(:,4)+NnOld,index(:,5)+NnOld+NedOld,index(:,3)+NnOld];

            % setting data of new mesh
            obj.setData;

        end

        %% tranforms and copy generations
        function mirror(obj, varargin)
            obj.nodes = ext_pmirror2(obj.nodes, varargin{:});
            obj.qcl = obj.qcl(:, [1, 4, 3, 2]);
            obj.clearSetDataFlag;
            obj.setData;
        end

        function newObj = getMirror(obj, varargin)
            newObj = copy(obj);
            newObj.nodes = ext_pmirror2(newObj.nodes, varargin{:});
            newObj.cl = newObj.cl(:, [1, 4, 3, 2]);
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

        function mzptr = getExtrude(obj, zLevels)
            % extrude quadrilateral mesh along z-axis and construct hexahedral mesh

            p2 = [obj.nodes, zeros(obj.Nn,1)];
            Np = size(p2,1);
            Nq = size(obj.qcl,1);
            Nz = numel(zLevels) - 1;

            if Nz < 1
                error('zLevels must contain at least two z coordinates');
            end

            % Build extruded node coordinates
            p3 = zeros(Np*(Nz+1), 3);
            for k = 1:(Nz+1)
                idx = (1:Np) + (k-1)*Np;
                p3(idx,1) = p2(:,1);
                p3(idx,2) = p2(:,2);
                p3(idx,3) = zLevels(k);
            end

            % Build hexahedral connectivity
            hhcl = zeros(Nq*Nz, 8);
            row = 1;

            for k = 1:Nz
                offsetBot = (k-1)*Np;
                offsetTop = k*Np;

                bot = obj.qcl + offsetBot;
                top = obj.qcl + offsetTop;

                % Reverse top face ordering for consistent hex orientation
                hhcl(row:row+Nq-1,:) = [bot(:,[1,2,3,4]), top(:,[1,2,3,4])];

                row = row + Nq;
            end

            mzptr = emdlab_m3d_hhmz(hhcl, p3);
            mzptr.color = obj.color;

        end

    end

    methods (Access = private)

        function clearSetDataFlag(obj)

            obj.isDataSet = false;
            obj.isAreaOfElementsEvaluated = false;
            obj.isMeshZoneAreaEvaluated = false;

        end

    end

end

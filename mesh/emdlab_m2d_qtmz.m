% EMDLAB: Electrical Machines Design Laboratory
% Quadrilateral-Triangular mesh zone (2D element)

classdef emdlab_m2d_qtmz <  handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_m2d_xmz

    methods
        %% Constructor and Destructor
        function obj = emdlab_m2d_qtmz(varargin)

            switch nargin
                case 2
                    obj.cl = varargin{1};
                    obj.nodes = varargin{2};

                case 3
                    obj.cl = [varargin{1}; varargin{2}, zeros(size(varargin{2},1),1)];
                    obj.nodes = varargin{3};

                otherwise
                    error('Wrong number of input arguments.');
            end

            obj = obj.setData;

        end

        %% Topological Functions
        function obj = setData(obj)

            % check if already data is set
            if obj.isDataSet, return; end

            % getting the number of elements
            ne = obj.Ne;

            % sort rows base on last column to distinquish
            % between quadrilateral and triangles
            obj.cl = sortrows(obj.cl,4);

            % number of triangles in mesh zone
            nt = sum(obj.cl(:,4) == 0);

            % number of quadrilaterals in mesh zone
            nq = ne - nt;

            % indices referring to triangles and quadrilaterals
            idx_t = 1:nt;
            idx_q = (nt+1):ne;

            % all edges
            e1t = obj.cl(idx_t,[1,2]);
            e2t = obj.cl(idx_t,[2,3]);
            e3t = obj.cl(idx_t,[3,1]);
            e1q = obj.cl(idx_q,[1,2]);
            e2q = obj.cl(idx_q,[2,3]);
            e3q = obj.cl(idx_q,[3,4]);
            e4q = obj.cl(idx_q,[4,1]);

            % stack all edges
            allEdges = [e1t; e2t; e3t; e1q; e2q; e3q; e4q];
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
            obj.elements = [reshape(ic(1:3*nt), nt, 3), zeros(nt,1); reshape(ic(3*nt+1:end), nq, 4)];

            % find boundary edges
            idx = emdlab_mex_findSignedPairs(obj.elements, size(obj.edges,1));
            obj.bedges = idx ~= 3;

            % change states
            obj.isDataSet = true;

        end

        %% Mesh Visiualization
        function varargout = showm(obj)

            f = emdlab_r2d_mesh();
            ax = axes(f);

            fcl = obj.cl(:,1:4);
            fcl(fcl == 0) = nan;

            patch('Faces',fcl,'Vertices',obj.nodes,'FaceColor',...
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

end

% EMDLAB: Electrical Machines Design Laboratory
% Quadrilateral mesh zone (2D element)

classdef emdlab_m2d_qmz <  handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_m2d_xmz

    methods
        %% Constructor and Destructor
        function obj = emdlab_m2d_qmz(cl, nodes)

            if nargin < 2, error('Not enough input arguments.'); end
            if nargin > 2, error('Too many input arguments.'); end

            obj.nodes = nodes;
            obj.cl = cl;
            obj = obj.setData;

        end

        %% Topological Functions        
        function obj = setData(obj)

            % check if already data is set
            if obj.isDataSet, return; end

            % getting the number of elements
            ne = obj.Ne;

            % quadrilateral edges
            e1 = obj.cl(:,[1,2]);
            e2 = obj.cl(:,[2,3]);
            e3 = obj.cl(:,[3,4]);
            e4 = obj.cl(:,[4,1]);

            % stack all edges
            allEdges = [e1; e2; e3; e4];
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
            obj.elements = reshape(ic, ne, 4);

            % find boundary edges
            idx = emdlab_mex_findSignedPairs(obj.elements, size(obj.edges,1));
            obj.bedges = idx ~= 3;

            % change states
            obj.isDataSet = true;

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
            NeOld = obj.Ne;
            NedOld = size(obj.edges, 1);
            % nodes of new mesh
            obj.nodes = [obj.nodes;...
                (obj.nodes(obj.edges(:,1),:)+obj.nodes(obj.edges(:,2),:))/2;
                (obj.nodes(obj.cl(:,1),:)+obj.nodes(obj.cl(:,2),:)+obj.nodes(obj.cl(:,3),:)+obj.nodes(obj.cl(:,4),:))/4];
            % index of nodes on old edges
            index = [abs(obj.elements),(1:NeOld)'];
            % new connctivity list
            obj.cl = [obj.cl(:,1),index(:,1)+NnOld,index(:,5)+NnOld+NedOld,index(:,4)+NnOld
                obj.cl(:,2),index(:,2)+NnOld,index(:,5)+NnOld+NedOld,index(:,1)+NnOld
                obj.cl(:,3),index(:,3)+NnOld,index(:,5)+NnOld+NedOld,index(:,2)+NnOld
                obj.cl(:,4),index(:,4)+NnOld,index(:,5)+NnOld+NedOld,index(:,3)+NnOld];

            % setting data of new mesh
            obj.setData;

        end

        

        

        function mzptr = getExtrude(obj, zLevels)
            % extrude quadrilateral mesh along z-axis and construct hexahedral mesh

            p2 = [obj.nodes, zeros(obj.Nn,1)];
            Np = size(p2,1);
            Nq = size(obj.cl,1);
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

                bot = obj.cl + offsetBot;
                top = obj.cl + offsetTop;

                % Reverse top face ordering for consistent hex orientation
                hhcl(row:row+Nq-1,:) = [bot(:,[1,2,3,4]), top(:,[1,2,3,4])];

                row = row + Nq;
            end

            mzptr = emdlab_m3d_hhmz(hhcl, p3);
            mzptr.color = obj.color;

        end

    end

end

% EMDLAB: Electrical Machines Design Laboratory
% Tetrahedral mesh zone (3D element)

classdef emdlab_m3d_thmz < handle & emdlab_g2d_constants & matlab.mixin.Copyable

    properties(SetAccess = private)

        % mesh nodes
        nodes (:,3) double;

        % mesh connectivity list
        cl (:,4) double;

        % mesh elements: [facet1, facet2, facet3, facet4]
        elements (:,4) double;

        % unique facets: [node1, node2, node3]
        facets (:,3) double;

        % list of boundary facets
        bfacets (:,1);

        % unique edges: [node1, node2]
        edges (:,2);

        % list of boundary edges
        bedges (:,1);

        % Named Selections
        nodeNamedSelections (1,1) struct;
        edgeNamedSelections (1,1) struct;
        facetNamedSelections (1,1) struct;

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

    end

     properties (Access = private)

        % element volume
        ev

        % mesh zone volume
        volume (1,1) double;
        
        % Q
        Q

        % Weight matrix
        Wm

        % states
        isDataSetted (1,1) logical = false;
        is_ev_Evaluated (1,1) logical = false;
        is_volume_Evaluated (1,1) logical = false;
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
        function obj = emdlab_m3d_thmz(cl, nodes)
            if nargin < 2, error('Not enough input arguments.'); end
            if nargin > 2, error('Too many input arguments.'); end
            obj.nodes = nodes;
            obj.cl = cl;
            obj.setdata;
        end

        function y = get.Nn(obj)
            y = size(obj.nodes,1);
        end

        function y = get.Ne(obj)
            y = size(obj.cl,1);
        end

        %% FEM preparation
        function evalVolumeOfElements(obj)
            if obj.is_ev_Evaluated, return; end
            % x, y and z coordinate of points
            xp = obj.nodes(:,1);
            yp = obj.nodes(:,2);
            zp = obj.nodes(:,3);
            % point coordinate of each triangle nodes
            xp = xp(obj.cl');
            yp = yp(obj.cl');
            zp = zp(obj.cl');
            % calculation the volume of each tetrahedras
            obj.ev = xp(1,:).*(yp(2,:).*(zp(4,:)-zp(3,:))+...
                yp(3,:).*(zp(2,:)-zp(4,:))+yp(4,:).*(zp(3,:)-zp(2,:))) +...
                xp(2,:).*(yp(1,:).*(zp(3,:)-zp(4,:))+...
                yp(3,:).*(zp(4,:)-zp(1,:))+yp(4,:).*(zp(1,:)-zp(3,:))) +...
                xp(3,:).*(yp(1,:).*(zp(4,:)-zp(2,:))+...
                yp(2,:).*(zp(1,:)-zp(4,:))+yp(4,:).*(zp(2,:)-zp(1,:))) +...
                xp(4,:).*(yp(1,:).*(zp(2,:)-zp(3,:))+...
                yp(2,:).*(zp(3,:)-zp(1,:))+yp(3,:).*(zp(1,:)-zp(2,:)));
            obj.ev = abs(obj.ev)/6;
            % change states
            obj.is_ev_Evaluated = true;
        end

        function evalVolume(obj)
            if obj.is_volume_Evaluated, return; end
            obj.evalVolumeOfElements;
            obj.volume = sum(obj.ev);
            % change states
            obj.is_volume_Evaluated = true;
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

            % first facet of each triangle
            f1 = obj.cl(:,[1,2,3]);
            % second facet of each triangle
            f2 = obj.cl(:,[2,4,3]);
            % third facet of each triangle
            f3 = obj.cl(:,[3,4,1]);
            % forth facet of each triangle
            f4 = obj.cl(:,[1,4,2]);

            % sorting for lower index
            [f1,s1] = sort(f1,2);
            [f2,s2] = sort(f2,2);
            [f3,s3] = sort(f3,2);
            [f4,s4] = sort(f4,2);

            % specefying changed facet index
            s1 = ((s1(:,1)==1)&(s1(:,2)==3))|...
                ((s1(:,1)==3)&(s1(:,2)==2))|...
                ((s1(:,1)==2)&(s1(:,2)==1));
            s2 = ((s2(:,1)==1)&(s2(:,2)==3))|...
                ((s2(:,1)==3)&(s2(:,2)==2))|...
                ((s2(:,1)==2)&(s2(:,2)==1));
            s3 = ((s3(:,1)==1)&(s3(:,2)==3))|...
                ((s3(:,1)==3)&(s3(:,2)==2))|...
                ((s3(:,1)==2)&(s3(:,2)==1));
            s4 = ((s4(:,1)==1)&(s4(:,2)==3))|...
                ((s4(:,1)==3)&(s4(:,2)==2))|...
                ((s4(:,1)==2)&(s4(:,2)==1));
            
            % unification of facets
            [obj.facets,~,ic] = unique([f1;f2;f3;f4],'rows');

            % getting number of elements
            ne = obj.Ne;

            % getting index of facets corresponding to each elements
            f1 = ic(1:ne);
            f2 = ic(1+ne:2*ne);
            f3 = ic(1+2*ne:3*ne);
            f4 = ic(1+3*ne:4*ne);

            % specefying boundary facets
            obj.bfacets = sparse([f1,f2,f3,f4],ones(4*ne,1),ones(4*ne,1));
            obj.bfacets = full(obj.bfacets == 1);
            
            % specefying trace direction
            f1(s1) = -f1(s1);
            f2(s2) = -f2(s2);
            f3(s3) = -f3(s3);
            f4(s4) = -f4(s4);
            
            % element matrix
            obj.elements = [f1,f2,f3,f4];
            
            % evaluation of area of each elements
            obj.evalVolumeOfElements;
            obj.facetNamedSelections.('none') = find(obj.bfacets);
            
            % change states
            obj.isDataSetted = true;

        end
       
        %% visiualization functions
        function showm(obj)
            obj.setdata;
            [f,ax] = emdlab_r3d_mesh();

            f.Name = ['[Global Mesh][','Nn = ',num2str(obj.Nn),'][Ne = ',num2str(obj.Ne),']'];
            patch('Faces',obj.facets(obj.bfacets,1:3),'Vertices',...
                obj.nodes,'FaceColor',...
                obj.color,'EdgeColor','k', 'parent', ax);
            set(f, 'Visible', 'on');
        end

        function showg(obj)
            obj.setdata;

            f = GraphicWindow;
            h = guihandles(f);

            patch('Faces',obj.facets(obj.bfacets,1:3),'Vertices',obj.nodes,...
                'FaceColor',obj.color,'EdgeColor','none','parent',h.va);

            tr = triangulation(obj.facets(obj.bfacets,1:3),obj.nodes);

            n = tr.faceNormal;
            ea = tr.edgeAttachments(tr.edges);
            e = [];

            for i = 1:numel(ea)
                tmp = sum(n(ea{i}(1),:).*n(ea{i}(2),:));
                if tmp>pi/180
                    e(end+1) = i;
                end
            end
            %       e = tr.featureEdges(pi/180);
            %

            eee = tr.edges;
            patch('Faces',eee(e,:),'Vertices',tr.Points,...
                'FaceColor','k','EdgeColor','k','parent',ah);
            set(gcf,'HandleVisibility','off', 'Visible', 'on');
        end

        function showwf(obj, color)
            if nargin<2
                color = 'c';
            end
            ah = setFigure(['[Global Mesh][','Nn = ',num2str(obj.Nn),'][Ne = ',num2str(obj.Ne),']'], true);
            patch('Faces',obj.facets(obj.bfacets,1:3),'Vertices',...
                obj.nodes,'FaceColor',...
                color,'EdgeColor','w',...
                'FaceAlpha',0.5, 'parent', ah);
            set(gcf,'HandleVisibility','off', 'Visible', 'on');
        end

        function shownf(obj, name)
            name = obj.checkFacetNamedSelectionExistence(name);
            axis off equal;
            patch('Faces',obj.facets(obj.facetNamedSelections.(name),1:3),'Vertices',...
                obj.nodes,'FaceColor',...
                'c','EdgeColor','b',...
                'FaceAlpha',1);
            set(gca,'Clipping','off');
            setFigure;
        end

        function shownfs(obj)
            axis off equal;
            hold all;
            nfs = fieldnames(obj.facetNamedSelections);
            for i = 1:numel(nfs)
                tmp = rand(1, 3);
                patch('Faces',obj.facets(obj.facetNamedSelections.(nfs{i}),1:3),'Vertices',...
                    obj.nodes,'FaceColor',...
                    tmp,'EdgeColor','k',...
                    'FaceAlpha',1);
            end
            legend(nfs);
            set(gca,'Clipping','off');
            setFigure;
        end
        
        %% Named Selections
        % node
        function name = checkNodeNamedSelectionExistence(obj,name)
            name = rmspaces(name);
            if ~isfield(obj.nodeNamedSelections,name)
                error('Specified node named selection does not exist.');
            end
        end
        function name = checkNodeNamedSelectionNonExistence(obj,name)
            name = rmspaces(name);
            if isfield(obj.nodeNamedSelections,name)
                error('Specified node named selection already exist.');
            end
        end
        function addNodeNamedSelection(obj, name, indices)
            name = obj.checkNodeNamedSelectionNonExistence(name);
            obj.nodeNamedSelections.(name) = indices;
        end
        % facet
        function name = checkFacetNamedSelectionExistence(obj,name)
            name = rmspaces(name);
            if ~isfield(obj.facetNamedSelections,name)
                error('Specified facet named selection does not exist.');
            end
        end
        function name = checkFacetNamedSelectionNonExistence(obj,name)
            name = rmspaces(name);
            if isfield(obj.facetNamedSelections,name)
                error('Specified facet named selection already exist.');
            end
        end
        function addFacetNamedSelection(obj, name, indices)
            name = obj.checkFacetNamedSelectionNonExistence(name);
            obj.facetNamedSelections.(name) = indices;
            obj.facetNamedSelections.('none') = setdiff(...
                obj.facetNamedSelections.('none'),...
                obj.facetNamedSelections.(name));
        end
        %% Tools Functions
        function moveNodes(obj,MovTol)
            if nargin<2
                MovTol = 1e-3;
            end
            % connectivity matrix for nodes
            Con = sparse(obj.edges(:,1),obj.edges(:,2),...
                ones(size(obj.edges,1),1),obj.Nn,obj.Nn);
            Con = Con + Con';
            % loop for movments
            inodes = obj.getinodes;
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
        %% Index Finding
        function y = getfbf(obj)
            y = find(obj.bfacets);
        end
        function y = getfbn(obj)
            y = obj.getfbf;
            y = unique(y(:));
        end
        function y = getbnodes(obj)
            % getting index of boundary nodes
            y = obj.facets(obj.bfacets,:);
            y = unique(y(:));
        end
        function y = getinodes(obj)
            % getting index of inner nodes
            y = obj.getbnodes;
            y = setdiff((1:obj.Nn)',y);
        end
        function y = getNodeIndexOnPlane(obj, p0, n)
            y = find(abs(obj.nodes*n'-p0*n') < obj.gleps);
        end
        function y = getNodeIndexOnHalfPlane(obj, p0, p1, p2)
            y = obj.getNodeIndexOnPlane(p0, cross(p1, p2));
            tmp = obj.nodes(y,:)*p1' >= 0;
            y = y(tmp);
        end
        function y = getFacetIndexOnPlane(obj, varargin)
            y = obj.getNodeIndexOnPlane(varargin{:});
            y = ismember(obj.facets(:,1),y)&...
                ismember(obj.facets(:,2),y)&ismember(obj.facets(:,3),y);
            y = find(y);
        end
        function y = getFacetIndexOnHalfPlane(obj, varargin)
            y = obj.getNodeIndexOnHalfPlane(varargin{:});
            y = ismember(obj.facets(:,1),y)&...
                ismember(obj.facets(:,2),y)&ismember(obj.facets(:,3),y);
            y = find(y);
        end

        %% geometrical operations
        function y = getCenterOfElements(obj)
            % get center of elements
            y = (obj.nodes(obj.cl(:,1),:) + ...
                obj.nodes(obj.cl(:,2),:) + ...
                obj.nodes(obj.cl(:,3),:)+...
                obj.nodes(obj.cl(:,4),:))/4;
        end
        
        function y = getVolumeOfElements(obj)
            obj.setdata;
            obj.evalVolumeOfElements;
            y = obj.ev;
        end
        
        function y = getVolume(obj)
            obj.setdata;
            obj.evalVolumeOfElements;
            obj.evalVolume;
            y = obj.volume;
        end
        
        function y = getSurfaceArea(obj)
            p1p2 = obj.nodes(obj.facets(obj.bfacets,2),:) - obj.nodes(obj.facets(obj.bfacets,1),:);
            p1p3 = obj.nodes(obj.facets(obj.bfacets,3),:) - obj.nodes(obj.facets(obj.bfacets,1),:);
            y = cross(p1p2, p1p3);
            y = 0.5*sum(sqrt(sum(y.^2,2)));
        end
        
        %% tranforms and copy generations
        function mirror(obj, p0, p1)

            % set default plane if p0 and p1 are not provided
            if nargin < 3
                p1 = p0;          % if only points provided, use p0 as reference for normal
                p0 = [0,0,0];     % default plane passes through origin
            end

            % input validation
            if ~isnumeric(p0) || ~isnumeric(p1)
                error('<p0> and <p1> must be numeric data.');
            elseif ~isequal(size(p0), [1,3]) || ~isequal(size(p1), [1,3])
                error('<p0> and <p1> must be 1x3 vectors.');
            end

            % plane normal vector
            n = p1 - p0;
            n = n / norm(n);
            obj.nodes = obj.nodes - 2 * ((obj.nodes - p0) * n') * n;

        end

        function newObj = getMirror(obj, varargin)
            newObj = copy(obj);
            newObj.mirror(varargin{:});
            newObj.cl = newObj.cl(:,[1,3,2,4]);
            newObj.makeFalse_isDataSetted;
            newObj.setdata;
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
            %ROTATE Rotate mesh nodes around an axis defined by (p0,p1) by angle rotAngle (rad)

            % set default axis if not provided
            if nargin < 4
                p1 = p0;
                p0 = [0,0,0];     % axis passes through origin if only one point given
            end

            % input validation
            if ~isnumeric(p0) || ~isnumeric(p1)
                error('<p0> and <p1> must be numeric data.');
            elseif ~isequal(size(p0), [1,3]) || ~isequal(size(p1), [1,3])
                error('<p0> and <p1> must be 1x3 vectors.');
            end

            % ---- Compute rotation axis ----
            u = p1 - p0;
            if norm(u) < eps
                error('Rotation axis is not defined. p0 and p1 must be distinct points.');
            end
            u = u / norm(u);

            % ---- Translate geometry ----
            V = obj.nodes - p0;   % shift nodes so axis goes through origin

            % ---- Rodrigues’ rotation ----
            cosT = cos(rotAngle);
            sinT = sin(rotAngle);

            dotUV   = V*u';                                      % Nx1
            crossUV = cross(repmat(u, size(V,1), 1), V, 2);      % Nx3

            V_rot = V*cosT + crossUV*sinT + dotUV*(1-cosT).*u;

            % ---- Translate back ----
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
            %SHIFT Shift mesh nodes so that p0 is mapped to p1

            % set default plane if p0 and p1 are not provided
            if nargin < 3
                p1 = p0;          % if only points provided, use p0 as reference for normal
                p0 = [0,0,0];     % default plane passes through origin
            end
            
            % input validation
            if ~isnumeric(p0) || ~isnumeric(p1)
                error('<p0> and <p1> must be numeric data.');
            elseif ~isequal(size(p0), [1,3]) || ~isequal(size(p1), [1,3])
                error('<p0> and <p1> must be 1x3 vectors.');
            end

            % displacement vector
            d = p1 - p0;

            % shift all nodes
            obj.nodes = obj.nodes + d;
        end

        function newObj = getShift(obj,varargin)
            newObj = copy(obj);
            newObj.shift(varargin{:});
        end
    
    end

    methods (Access = private)

        function makeFalse_isDataSetted(obj)
            obj.isDataSetted = false;
            obj.is_ev_Evaluated = false;
            obj.is_volume_Evaluated = false;
            obj.is_Wm_Evaluated = false;
            obj.is_Q_Evaluated = false;
        end

    end

end

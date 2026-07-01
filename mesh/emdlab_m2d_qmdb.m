% EMDLAB: Electrical Machines Design Laboratory
% 2D quadrilateral mesh data base class

classdef emdlab_m2d_qmdb < handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_mdb_cp

    properties (SetAccess = private)

        % mesh nodes: [x,y]
        nodes (:,2) double;

        % mesh connectivity list
        cl (:,:) double;

        % mesh elements: [edge1, edge2, edge3, edge4, zone index]
        elements (:,5) double;

        % unique edges
        edges

        % list of boundary edges
        bedges

        % jacobian inverse transpose
        JIT

        % global element area
        gea (1,:) double;

        % element zone index
        ezi (:,:) logical;

        % elements material index
        emi (:,:) logical;

        % auxiliary stored matricies
        mtcs (1,1) struct;

        % named selections
        edgeNamedSelections (1,1) struct;

        % flag to print the elapsed times
        printFlag (1,1) logical = true;

        % element type
        etype (1,:) char = 'QL4';

    end

    properties (Dependent = true)

        % Number of nodes
        Nn

        % Number of elements
        Ne

        % flags for element type
        isQL4 (1,1) logical;
        isQL8 (1,1) logical;

    end

    properties (Access = private)

        % element type
        isd2ElementsGenerated (1,1) logical = false;
        isd3ElementsGenerated (1,1) logical = false;
        isJITEvaluated (1,1) logical = false;
        isKeFe_TL3_Evaluated (1,1) logical = false;
        isKexy4Fe_TL3_Evaluated (1,1) logical = false;
        isKeFe_TL6_Evaluated (1,1) logical = false;
        isKeFe_TL6_cte_Evaluated (1,1) logical = false;
        isKe_TL3_Evaluated (1,1) logical = false;
        isFe_TL3_Evaluated (1,1) logical = false;
        isMe_TL3_Evaluated (1,1) logical = false;
        isKexy1_TL3_Evaluated (1,1) logical = false;
        isKexy4_TL3_Evaluated (1,1) logical = false;
        isKexy1_TL6_Evaluated (1,1) logical = false;
        isKerz1_TL3_Evaluated (1,1) logical = false;

    end

    methods
        %% Constructor and Destructor
        function obj = emdlab_m2d_qmdb()
            % add default material
            obj.addMaterial('air');
        end

        function addMeshZone(obj, varargin)

            % you can pass an <emdlab_m2d_qmz> class without name
            if nargin == 2

                if ~isa(varargin{1}, 'emdlab_m2d_qmz')
                    error('Mesh zone class must be <emdlab_m2d_qmz>.');
                end
                mzName = obj.getDefaultMeshZoneName;
                mzptr = varargin{1};

                % you can pass an <emdlab_m2d_qmz> mesh zone with specified name
            elseif nargin == 3

                mzName = obj.checkMeshZoneNonExistence(varargin{1});
                if ~isa(varargin{2}, 'emdlab_m2d_qmz')
                    error('Mesh zone class must be <emdlab_m2d_qmz>.');
                end
                mzptr = varargin{2};

            else
                error('Wrong number of arguments.');
            end

            % adding new mesh zone
            obj.mzs.(mzName) = mzptr;
            obj.mzs.(mzName).material = 'air';
            obj.mzs.(mzName).color = rand(1,3);

            % changing states
            obj.makeFalse_isGlobalMeshGenerated;

        end

        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end

        function y = get.Ne(obj)
            y = size(obj.cl, 1);
        end

        function delete(obj)

            mzNames = fieldnames(obj.mzs);

            for i = 1:numel(mzNames)
                delete(obj.mzs.(mzNames{i}));
            end

        end

        %% FEM Preparation
        % generate global mesh
        function ggmesh(obj)

            % check states
            if obj.isGlobalMeshGenerated, return; end

            % generation of initial mesh
            Nn_tmp = 0;
            Ne_tmp = 0;
            mzNames = obj.getMeshZoneNames;

            for mzName = mzNames
                Nn_tmp = Nn_tmp + obj.mzs.(mzName).Nn;
                Ne_tmp = Ne_tmp + obj.mzs.(mzName).Ne;
            end

            % initialization of nodes and elements
            obj.nodes = zeros(Nn_tmp,2);
            obj.cl = zeros(Ne_tmp,4);
            obj.elements = zeros(Ne_tmp,5);
            nindex = 0;
            eindex= 0;

            % loop over mesh zones: for insertion of nodes and elements
            for i = 1:numel(mzNames)

                % insertion of nodes
                obj.nodes(1 + nindex:obj.mzs.(mzNames(i)).Nn + nindex, :) = ...
                    obj.mzs.(mzNames(i)).nodes;

                % insertion of elements
                obj.cl(1 + eindex:obj.mzs.(mzNames(i)).Ne + eindex, 1:4) = ...
                    obj.mzs.(mzNames(i)).cl + nindex;

                % specefying zone index
                obj.elements(1 + eindex:obj.mzs.(mzNames(i)).Ne + eindex, 5) = i;
                obj.mzs.(mzNames(i)).zi = i;
                nindex = nindex + obj.mzs.(mzNames(i)).Nn;
                eindex = eindex + obj.mzs.(mzNames(i)).Ne;

            end

            % unification of nodes
            [obj.nodes, ~, ic] = uniquetol(obj.nodes, obj.gleps, 'ByRows', true);
            obj.cl = ic(obj.cl);

            % setting l2g
            nindex = 0;

            % loop over mesh zones: for setting l2g
            for mz = mzNames
                obj.mzs.(mz).l2g = ic(nindex + 1:nindex + ...
                    obj.mzs.(mz).Nn);
                nindex = nindex + obj.mzs.(mz).Nn;        
            end

            obj.setData;
            obj.evalezi;

            % change states
            obj.isGlobalMeshGenerated = true;

            % settig element type
            obj.etype = 'QL4';

        end
        
        % evaluate element material & zone index
        function evalezi(obj)

            % construct emi & ezi
            meshZoneNames = obj.getMeshZoneNames;
            obj.materialNames = string(fieldnames(obj.mts)');
            obj.emi = false(obj.Nmts, obj.Ne);
            obj.ezi = false(obj.Ne, obj.Nmzs);

            for i = 1:obj.Nmzs
                obj.ezi(:,i) = obj.elements(:,5) == i;
                for j = 1:obj.Nmts
                    if obj.mzs.(meshZoneNames(i)).material == obj.materialNames(j)
                        obj.emi(j,obj.ezi(:,i)) = true;
                    end
                end
            end

        end
        
        function evalJIT(obj)
            % check states
            if obj.isJITEvaluated, return; end
            tic, disp('-------------------------------------------------------');
            % prerequisite
            obj.ggmesh;
            % getting number of elements
            xNe = obj.Ne;
            % forming gea vector
            obj.gea = zeros(1,xNe);
            mznames = fieldnames(obj.mzs);
            for i = 1:obj.Nmzs
                obj.gea(obj.ezi(:,obj.mzs.(mznames{i}).zi)) = ...
                    obj.mzs.(mznames{i}).getAreaOfElements;
            end
            % evaluation of jacobian inverse using d1 element data
            obj.JIT = tmdbc_evalJIT(obj.cl, obj.nodes, obj.gea);
            % change states
            obj.isJITEvaluated = true;
            disp('Evaluation of JIT completed.');
            toc, disp('-------------------------------------------------------');
        end
        
        function evalKeFe_TL3(obj)
            obj.ggmesh;
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            if obj.isKeFe_TL3_Evaluated, return; end
            % getting number of elements
            xNe = obj.Ne;
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy] = tmdbc_evalGphixy_TL3(obj.JIT);
            % evaluation of Ke
            obj.Ke = zeros(6,xNe);
            temp = 0;
            for i = 1:3
                for j = 1:i
                    temp = temp + 1;
                    obj.Ke(temp,:) = obj.gphix(i,:).*obj.gphix(j,:) +...
                        obj.gphiy(i,:).*obj.gphiy(j,:);
                end
            end
            % multiplying by triangle areas
            obj.Ke = obj.Ke*sparse(1:xNe,1:xNe,obj.gea);
            % elemental force matrix
            obj.Fe = repmat((obj.gea/3),3,1);
            disp('Calculation of [Ke] and [Fe] completed.');
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKeFe_TL3_Evaluated = true;
        end
        
        function evalKexy1_TL3(obj)
            if obj.isKexy1_TL3_Evaluated, return; end
            obj.ggmesh;
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            % getting number of elements
            xNe = obj.Ne;
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy] = tmdbc_evalGphixy_TL3(obj.JIT);
            % evaluation of Ke
            obj.Ke = zeros(6,xNe);
            temp = 0;
            for i = 1:3
                for j = 1:i
                    temp = temp + 1;
                    obj.Ke(temp,:) = (obj.gphix(i,:).*obj.gphix(j,:) +...
                        obj.gphiy(i,:).*obj.gphiy(j,:)).*obj.gea;
                end
            end
            disp('Calculation of [Kexy1] completed.');
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexy1_TL3_Evaluated = true;
        end
        
        function evalKexy4Fe_TL3(obj)
            obj.ggmesh;
            obj.evalJIT;
            if obj.isKexy4Fe_TL3_Evaluated
                return
            end
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % getting specefied edata
            edata = getedata('TL3');
            % evaluation of grad phi i
            gphi = obj.JIT*repmat(edata.GG,xNe,1);
            % inserting in gradx and grady phi i
            obj.gphix = (gphi(1:2:end,:))';
            obj.gphiy = (gphi(2:2:end,:))';
            % evaluation of Ke
            obj.Ke.x = zeros(6,xNe);
            obj.Ke.y = zeros(6,xNe);
            obj.Ke.xy = zeros(9,xNe);
            obj.Ke.yx = zeros(9,xNe);
            temp = 0;
            for i = 1:3
                for j = 1:i
                    temp = temp + 1;
                    obj.Ke.x(temp,:) = obj.gphix(i,:).*obj.gphix(j,:);
                    obj.Ke.y(temp,:) = obj.gphiy(i,:).*obj.gphiy(j,:);
                end
            end
            temp = 0;
            for i = 1:3
                for j = 1:3
                    temp = temp + 1;
                    obj.Ke.xy(temp,:) = obj.gphix(i,:).*obj.gphiy(j,:);
                    obj.Ke.yx(temp,:) = obj.gphiy(i,:).*obj.gphix(j,:);
                end
            end
            % multiplying by triangle areas
            obj.Ke.x = obj.Ke.x*sparse(1:xNe,1:xNe,obj.gea);
            obj.Ke.y = obj.Ke.y*sparse(1:xNe,1:xNe,obj.gea);
            obj.Ke.xy = obj.Ke.xy*sparse(1:xNe,1:xNe,obj.gea);
            obj.Ke.yx = obj.Ke.yx*sparse(1:xNe,1:xNe,obj.gea);
            % elemental force matrix
            obj.Fe = repmat((obj.gea/3),3,1);
            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
            % change states
            obj.isKexy4Fe_TL3_Evaluated = true;
        end
        
        function evalKeFe_TL6_cte(obj)
            obj.ggmesh;
            evalJIT;
            if obj.isKeFe_TL6_cte_Evaluated
                return
            end
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % getting specefied edata
            edata = getedata('TL6');
            % GG at three point of TR
            % at first point
            GG120 = repmat(edata.GG(0.5,0),xNe,1);
            % at second point
            GG1212 = repmat(edata.GG(0.5,0.5),xNe,1);
            % at third point
            GG012 = repmat(edata.GG(0,0.5),xNe,1);
            % evaluation of Ke
            obj.Ke = zeros(21,xNe);
            temp = 0;
            for i = 1:6
                for j = i:6
                    temp = temp + 1;
                    ke = (obj.JIT*GG120(:,j)).*(obj.JIT*GG120(:,i)) +...
                        (obj.JIT*GG1212(:,j)).*(obj.JIT*GG1212(:,i)) +...
                        (obj.JIT*GG012(:,j)).*(obj.JIT*GG012(:,i)) ;
                    ke = reshape(ke,2,[]);
                    obj.Ke(temp,:) = sum(ke);
                end
            end
            % multiplying by triangle areas
            obj.Ke = obj.Ke*sparse(1:xNe,1:xNe,obj.gea/3);
            % calculation of fi
            obj.Fe = [zeros(3,xNe);repmat((obj.gea/3),3,1)];
            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
            % change states
            obj.isKeFe_TL6_cte_Evaluated = true;
        end
        
        function evalKeFe_TL6(obj)
            obj.ggmesh;
            obj.evalJIT;
            if obj.isKeFe_TL6_Evaluated
                return
            end
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % getting specefied edata
            edata = getedata('TL6');
            % JITGG at six point of TR
            % point1
            tmp = obj.JIT*repmat(edata.GG(0,0),xNe,1);
            obj.Aux.JITGG00x = tmp(1:2:end,:)';
            obj.Aux.JITGG00y = tmp(2:2:end,:)';
            % point 2
            tmp = obj.JIT*repmat(edata.GG(1,0),xNe,1);
            obj.Aux.JITGG10x = tmp(1:2:end,:)';
            obj.Aux.JITGG10y = tmp(2:2:end,:)';
            % point 3
            tmp = obj.JIT*repmat(edata.GG(0,1),xNe,1);
            obj.Aux.JITGG01x = tmp(1:2:end,:)';
            obj.Aux.JITGG01y = tmp(2:2:end,:)';
            % point 4
            tmp = obj.JIT*repmat(edata.GG(0.5,0),xNe,1);
            obj.Aux.JITGG120x = tmp(1:2:end,:)';
            obj.Aux.JITGG120y = tmp(2:2:end,:)';
            % point 5
            tmp = obj.JIT*repmat(edata.GG(0.5,0.5),xNe,1);
            obj.Aux.JITGG1212x = tmp(1:2:end,:)';
            obj.Aux.JITGG1212y = tmp(2:2:end,:)';
            % point 6
            tmp = obj.JIT*repmat(edata.GG(0,0.5),xNe,1);
            obj.Aux.JITGG012x = tmp(1:2:end,:)';
            obj.Aux.JITGG012y = tmp(2:2:end,:)';
            % evaluation of Ke
            obj.Ke = zeros(21,xNe);
            tmp = 0;
            for i = 1:6
                for j = 1:i
                    tmp = tmp + 1;
                    obj.Ke(tmp,:) = (obj.Aux.JITGG120x(j,:)).*(obj.Aux.JITGG120x(i,:)) + ...
                        (obj.Aux.JITGG120y(j,:)).*(obj.Aux.JITGG120y(i,:)) +...
                        (obj.Aux.JITGG1212x(j,:)).*(obj.Aux.JITGG1212x(i,:)) + ...
                        (obj.Aux.JITGG1212y(j,:)).*(obj.Aux.JITGG1212y(i,:)) +...
                        (obj.Aux.JITGG012x(j,:)).*(obj.Aux.JITGG012x(i,:)) + ...
                        (obj.Aux.JITGG012y(j,:)).*(obj.Aux.JITGG012y(i,:));
                end
            end
            % multiplying by triangle areas
            obj.Ke = obj.Ke*sparse(1:xNe,1:xNe,obj.gea/3);
            % calculation of fi
            obj.Fe = [zeros(3,xNe);repmat((obj.gea/3),3,1)];
            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
            % change states
            obj.isKeFe_TL6_Evaluated = true;
        end
        
        function evalKeFe(obj,etype)
            obj.etype = etype;
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % evaluation of jacobian inverse using d1 element data
            edata = getedata('TL3');
            % x and y coordinate of points
            xp = obj.nodes(:,1);
            yp = obj.nodes(:,2);
            % point coordinate of each triangle nodes
            xp = xp(obj.cl(:,1:3)');
            yp = yp(obj.cl(:,1:3)');
            % calculation of area of each triangles
            obj.gea = zeros(1,xNe);
            mznames = fieldnames(obj.mzs);
            for i = 1:obj.Nmzs
                obj.gea(obj.ezi(:,obj.mzs.(mznames{i}).zi)) = ...
                    obj.mzs.(mznames{i}).getAreaOfElements;
            end
            % evaluation of a and b
            acoefs = edata.M\xp;
            bcoefs = edata.M\yp;
            % Jacobian inverse transpose of each elements
            obj.JIT = [bcoefs(3,:);-acoefs(3,:);...
                -bcoefs(2,:);acoefs(2,:)];
            obj.JIT = obj.JIT*sparse(1:xNe,1:xNe,1./(2*obj.gea));
            [i,j] = getij(2,xNe);
            obj.JIT = sparse(i,j,obj.JIT(:));
            % getting specefied edata
            edata = getedata(obj.etype);
            switch upper(obj.etype)
                case 'TL3'
                    % evaluation of grad phi i
                    gphi = obj.JIT*repmat(edata.GG,xNe,1);
                    % inserting in gradx and grady phi i
                    obj.gphix = (gphi(1:2:end,:))';
                    obj.gphiy = (gphi(2:2:end,:))';
                    % evaluation of Ke
                    obj.Ke = zeros(6,xNe);
                    temp = 0;
                    for i = 1:3
                        for j = 1:i
                            temp = temp + 1;
                            obj.Ke(temp,:) = obj.gphix(i,:).*obj.gphix(j,:) +...
                                obj.gphiy(i,:).*obj.gphiy(j,:);
                        end
                    end
                    % multiplying by triangle areas
                    obj.Ke = obj.Ke*sparse(1:xNe,1:xNe,obj.gea);
                    % elemental force matrix
                    obj.Fe = repmat((obj.gea/3),3,1);
                case 'TL6'
                    % GG at three point of TR
                    % at first point
                    GG120 = repmat(edata.GG(0.5,0),xNe,1);
                    % at second point
                    GG1212 = repmat(edata.GG(0.5,0.5),xNe,1);
                    % at third point
                    GG012 = repmat(edata.GG(0,0.5),xNe,1);
                    % evaluation of Ke
                    obj.Ke = zeros(21,xNe);
                    temp = 0;
                    for i = 1:6
                        for j = i:6
                            temp = temp + 1;
                            ke = (obj.JIT*GG120(:,j)).*(obj.JIT*GG120(:,i)) +...
                                (obj.JIT*GG1212(:,j)).*(obj.JIT*GG1212(:,i)) +...
                                (obj.JIT*GG012(:,j)).*(obj.JIT*GG012(:,i)) ;
                            ke = reshape(ke,2,[]);
                            obj.Ke(temp,:) = sum(ke);
                        end
                    end
                    % multiplying by triangle areas
                    obj.Ke = obj.Ke*sparse(1:xNe,1:xNe,obj.gea/3);
                    % calculation of fi
                    obj.Fe = [zeros(3,xNe);repmat((obj.gea/3),3,1)];
            end
            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
        end
        
        function evalMe_TL3(obj)
            if obj.isMe_TL3_Evaluated, return; end
            tic, disp('-------------------------------------------------------');
            obj.Me = repmat((1/6)*[1;0.5;1;0.5;0.5;1],1,obj.Ne)...
                *sparse(1:obj.Ne,1:obj.Ne,obj.gea);
            disp('Calculation of [Me] completed.');
            toc, disp('-------------------------------------------------------');
            % changing states
            obj.isMe_TL3_Evaluated = true;
        end

        %% Higher Order Elements
        function gd2elements(obj)
            if obj.isd2ElementsGenerated, return; end
            obj.ggmesh;
            if obj.etype ~= 'TL3'
                obj.makeFalse_isGlobalMeshGenerated;
                obj.ggmesh;
            end
            Nn_old = obj.Nn;
            obj.nodes = [obj.nodes;...
                (obj.nodes(obj.edges(:,1),:)+obj.nodes(obj.edges(:,2),:))/2];
            obj.cl = [obj.cl,abs(obj.elements(:,1:3))+Nn_old];
            % settig element type
            obj.etype = 'TL6';
            % change states
            obj.isd2ElementsGenerated = true;
        end
        
        function gd3elements(obj)
            % uncompleted
            if obj.isd3ElementsGenerated
                return
            end
            obj.ggmesh;
            if obj.etype ~= 'TL3'
                obj.makeFalse_isGlobalMeshGenerated;
                obj.ggmesh;
            end
            Nn_old = obj.Nn;
            obj.nodes = [obj.nodes;...
                obj.nodes(obj.edges(:,1),:)*2/3+obj.nodes(obj.edges(:,2),:)/3;...
                ];
            obj.cl = [obj.cl,abs(obj.elements(:,1:3))+Nn_old];
            % settig element type
            obj.etype = 'TL6';
            % change states
            obj.isd3ElementsGenerated = true;
        end
        
        %% Topological Functions
        % setting needed data
        function setData(obj)

            % quadrilateral edges
            e1 = obj.cl(:,[1,2]);
            e2 = obj.cl(:,[2,3]);
            e3 = obj.cl(:,[3,4]);
            e4 = obj.cl(:,[4,1]);

            % sorting for lower index
            [e1,s1] = sort(e1,2);
            [e2,s2] = sort(e2,2);
            [e3,s3] = sort(e3,2);
            [e4,s4] = sort(e4,2);

            % specefying changed edge index
            s1 = s1(:, 1) == 2;
            s2 = s2(:, 1) == 2;
            s3 = s3(:, 1) == 2;
            s4 = s4(:, 1) == 2;

            % unification of edges
            [obj.edges, ~, ic] = unique([e1; e2; e3; e4], 'rows');
            
            % getting number of elements
            ne = obj.Ne;

            % getting index of edge corresponding to each elements
            e1 = ic(1:ne);
            e2 = ic(1 + ne:2 * ne);
            e3 = ic(1 + 2 * ne:3 * ne);
            e4 = ic(1 + 3 * ne:4 * ne);

             % specefying boundary edges
            obj.bedges = sparse([e1, e2, e3, e4], ones(4 * ne, 1), ones(4 * ne, 1));
            obj.bedges = full(obj.bedges == 1);

            % specefying trace direction
            e1(s1) = -e1(s1);
            e2(s2) = -e2(s2);
            e3(s3) = -e3(s3);
            e4(s4) = -e4(s4);

            % element matrix
            obj.elements(:, 1:4) = [e1, e2, e3, e4];

            % edge element
            obj.edges = [obj.edges, zeros(size(obj.edges,1), 6)];
            emdlab_m2d_qmdbc_evalee(obj.edges, obj.elements);

        end

        function strefine(obj)
            mznames = fieldnames(obj.mzs);
            for i = 1:numel(mznames)
                obj.mzs.(mznames{i}).strefine;
            end
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function adrefine(obj, ti)
            el = sqrt(sum((obj.nodes(obj.edges(:,1),:) - ...
                obj.nodes(obj.edges(:,2),:)).^2,2));
            b1 = el(abs(obj.elements(:,1)));
            b2 = el(abs(obj.elements(:,2)));
            b3 = el(abs(obj.elements(:,3)));
            [~,baseEdge] = max([b1, b2, b3], [], 2);
            for i = 1:length(ti)
                bisectTriangle(ti(i));
            end
            function bisectTriangle(index)
                eIndex = abs(obj.elements(index, baseEdge(index)));
                s = sign(obj.elements(index,:));
                clRow = obj.cl(index, :);
                eRow = obj.edges(eIndex, :);
                if obj.bedges(eIndex)
                    ptmp = (obj.nodes(eRow(1),:)+obj.nodes(eRow(2),:))/2;
                    obj.nodes(end+1,:) = ptmp;
                    nIndex = obj.Nn;
                    switch baseEdge(index)
                        case 1
                            %                   clRow1 =  [clRow(1), nIndex]
                        case 2
                        case 3
                    end

                    e1 = eRow;
                    e1(2) = nIndex;
                    e2 = eRow;
                    e2(1) = nIndex;
                    e3 = eRow;
                    disp('salam');
                else

                end
            end
        end
        %% mesh visiualization
        % show global mesh
        function varargout = showm(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:}); 
            if isa(f,"matlab.ui.Figure")
                f.MenuBar = "none";
            end
            
            obj.ggmesh;
            mzNames = string(fieldnames(obj.mzs)');

            for mzName = mzNames
                plt = patch(ax,'Faces', obj.mzs.(mzName).cl, ...
                    'Vertices', obj.mzs.(mzName).nodes, 'FaceColor', ...
                    'c', 'EdgeColor', [0.2, 0.2, 0.2], ...
                    'FaceAlpha', 0.7, ...
                    'HitTest','on','PickableParts','visible');
                plt.UserData = mzName;
            end

            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch(ax,'Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5,'HitTest','off','PickableParts','none');

            zoom on;box on;
            ax.Color = [0.86,0.86,0.86];
            grid on;
            grid minor;
            axis(ax, 'equal');
            ax.Toolbar.Visible = 'off';
            ax.GridColor      = [0.7 0.7 0.7];
            ax.MinorGridColor = [0.5 0.5 0.5];
            ax.GridAlpha      = 1;
            ax.MinorGridAlpha = 1;

            set(gcf,'WindowButtonMotionFcn',@hoverFcn);

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

            function hoverFcn(src,~)
                h = hittest(src);
                for i = 1:numel(ax.Children)
                    if isequal(h,ax.Children(i))
                        if isa(ax.Children(i), 'matlab.graphics.primitive.Patch')
                            e.Button = 1;
                            emdlab_flib_selectPatchCallbackGM(ax.Children(i),e);
                            return;
                        end
                    end
                end
                for i = 1:numel(ax.Children)
                    if isa(ax.Children(i), 'matlab.graphics.primitive.Patch')
                        if ischar(ax.Children(i).FaceColor)
                            if strcmpi(ax.Children(i).FaceColor, 'c')
                                set(ax.Children(i), 'FaceColor', 'c', 'FaceAlpha', 0.7);
                                drawnow;
                            end
                        else
                            if any(ax.Children(i).FaceColor ~= [0,1,1])
                                set(ax.Children(i), 'FaceColor', 'c', 'FaceAlpha', 0.7);
                                drawnow;
                            end
                        end
                    end
                end
                title(ax,'');
            end

        end

        function varargout = showgg(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:}); 
            if isa(f,"matlab.ui.Figure")
                f.MenuBar = "none";
            end
            
            obj.ggmesh;
            mzNames = string(fieldnames(obj.mzs)');

            for mzName = mzNames
                plt = patch(ax,'Faces', obj.mzs.(mzName).cl, ...
                    'Vertices', obj.mzs.(mzName).nodes, 'FaceColor', ...
                    'c', 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'HitTest','on','PickableParts','visible');
                plt.UserData = mzName;
            end

            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch(ax,'Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5,'HitTest','off','PickableParts','none');

            zoom on;box on;
            ax.Color = [0.86,0.86,0.86];
            grid on;
            grid minor;
            axis(ax, 'equal');
            ax.Toolbar.Visible = 'off';
            ax.GridColor      = [0.4 0.4 0.4];
            ax.MinorGridColor = [0.2 0.2 0.2];
            ax.GridAlpha      = 1;
            ax.MinorGridAlpha = 1;

            set(gcf,'WindowButtonMotionFcn',@hoverFcn);

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

            function hoverFcn(src,~)
                h = hittest(src);
                for i = 1:numel(ax.Children)
                    if isequal(h,ax.Children(i))
                        if isa(ax.Children(i), 'matlab.graphics.primitive.Patch')
                            e.Button = 1;
                            emdlab_flib_selectPatchCallbackGM(ax.Children(i),e);
                            updateXLabel;
                            return;
                        end
                    end
                end
                for i = 1:numel(ax.Children)
                    if isa(ax.Children(i), 'matlab.graphics.primitive.Patch')
                        if ischar(ax.Children(i).FaceColor)
                            if strcmpi(ax.Children(i).FaceColor, 'c')
                                set(ax.Children(i), 'FaceColor', 'c', 'FaceAlpha', 0.5);
                                drawnow;
                            end
                        else
                            if any(ax.Children(i).FaceColor ~= [0,1,1])
                                set(ax.Children(i), 'FaceColor', 'c', 'FaceAlpha', 0.5);
                                drawnow;
                            end
                        end
                    end
                end
                title(ax,'');
                updateXLabel;

                function updateXLabel()
                    % Get cursor position in axes units
                    cp = ax.CurrentPoint;
                    x = cp(1,1);
                    y = cp(1,2);

                    % Check if cursor is inside axes limits
                    xl = ax.XLim;
                    yl = ax.YLim;

                    if x < xl(1) || x > xl(2) || y < yl(1) || y > yl(2)
                        return
                    end

                    % Update xlabel
                    xlabel(ax, sprintf('X = %.2f ,  Y = %.2f ,  R = %.2f ,  D = %.2f',...
                        x, y, norm([x,y]), 2*norm([x,y])), 'Interpreter','none');
                end

            end

        end
        
        % show free boundary
        function varargout = showfb(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;
%             title(ax,['Global mesh free boundary edges', ', Nbe = ', num2str(sum(obj.bedges))]);

            patch('Faces', obj.edges(obj.bedges, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end
        
        % show wire frame
        function varargout = showwf(obj, varargin)

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

        % show mesh zones
        function varargout = showmzs(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            mzNames = fieldnames(obj.mzs);
            %             title(ax,['Number of mesh zones = ', num2str(numel(mzNames))]);

            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                patch('Faces', mzptr.cl, ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    mzptr.color, 'linewidth', 0.05 ,'EdgeColor', [0, 0, 0], ...
                    'FaceAlpha', 1, 'Parent', ax);
            end

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function showmz(obj, mzname)
            mzname = obj.checkMeshZoneExistence(mzname);
            obj.showm('c');
            hold on;
            patch('Faces',obj.mzs.(mzname).cl,...
                'Vertices',obj.mzs.(mzname).nodes,'FaceColor','r','EdgeColor','r',...
                'FaceAlpha',0.1);
        end
        
        function showmd(obj, color)
            if ~obj.isGlobalMeshGenerated
                error('Global mesh is not generated.');
            end
            if nargin<2
                color = 'k';
            end
            ah =  setFigure(['[Global Mesh][','Nn = ',num2str(obj.Nn),'][Ne = ',num2str(obj.Ne),']']);
            patch('Faces',obj.cl(:,1:3),'Vertices',obj.nodes,'FaceColor',...
                'c','FaceAlpha',0.5,'EdgeColor',color,'parent',ah);
            plot(obj.nodes(:,1),obj.nodes(:,2),...
                'o','color','k','linewidth',1,'MarkerFaceColor','k', 'parent', ah);
            axis(ah,'off','equal')
            set(gcf,'HandleVisibility','off');
        end

        %% tools: copy and transform
        % copy and transform
        function copyMirrorMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getMirror(varargin{:});
        end

        function cmmz(varargin)
            copyMirrorMeshZone(varargin{:});
        end

        function copyRotateMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getRotate(varargin{:});
        end

        function crmz(varargin)
            copyRotateMeshZone(varargin{:});
        end

        function copyShiftMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getShift(varargin{:});
        end

        function cshmz(varargin)
            copyShiftMeshZone(varargin{:});
        end

        % only transform
        function mirrorMeshZone(obj, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(mzName).mirror(varargin{:});
        end

        function mmz(varargin)
            mirrorMeshZone(varargin{:});
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

        %% Tools: some operations on mesh zones
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
        function getQuality(obj)
            % edges length
            el = sqrt(sum((obj.nodes(obj.edges(:,1),:) - ...
                obj.nodes(obj.edges(:,2),:)).^2,2));
            b1 = el(abs(obj.elements(:,1)));
            b2 = el(abs(obj.elements(:,2)));
            b3 = el(abs(obj.elements(:,3)));
            % mesh quality
            y = ((b1+b2-b3).*(b1-b2+b3).*(-b1+b2+b3))./(b1.*b2.*b3);
            fprintf('Average Quality = %f\n',mean(y));
            fprintf('Minimum Quality = %f\n',min(y));
        end
        function gq(varargin)
            getQuality(varargin{:});
        end
        function ttmdbc = getExtrude(obj, z, skewAngle)
            ttmdbc = TTMDBC();
            mzNames = fieldnames(obj.mzs);
            if iscolumn(z)
                z = z';
            end
            Nz = length(z);
            if nargin<3
                for i = 1:numel(mzNames)
                    mzptr = obj.mzs.(mzNames{i});
                    ztmp = repmat(z, mzptr.Nn, 1);
                    ttmdbc.addmz(mzNames{i}, TTMZPC(tmzpc_getExtrude(...
                        mzptr.cl,obj.elements(obj.ezi(:,mzptr.zi),1:3),...
                        mzptr.Nn,Nz-1),[repmat(mzptr.nodes,Nz,1),ztmp(:)]));
                    ttmdbc.mzs.(mzNames{i}).material = mzptr.material;
                end
            else
                stepAngle = skewAngle/(Nz-1);
                p = zeros(obj.Nn*Nz,3);
                for i = 1:Nz
                    p((i-1)*obj.Nn+1:i*obj.Nn,1:3) = ...
                        ext_protate3z([obj.nodes,repmat(z(i),obj.Nn,1)],(i-1)*stepAngle);
                end
                ttmz = TTMZPC(tmzpc_getExtrude(obj.cl,obj.elements,...
                    obj.Nn,Nz-1),p);
            end
        end
        function getMeshZoneExtrude(obj, ttmptr, mzName, z, skewAngle)
            mzName = obj.checkMeshZoneExistence(mzName);
            if iscolumn(z)
                z = z';
            end
            Nz = length(z);
            mzptr = obj.mzs.(mzName);
            ztmp = repmat(z, mzptr.Nn, 1);
            if nargin<5
                ttmptr.addmz(mzName, TTMZPC(tmzpc_getExtrude(...
                    mzptr.cl,obj.elements(obj.ezi(:,mzptr.zi),1:3),...
                    mzptr.Nn,Nz-1),[repmat(mzptr.nodes,Nz,1),ztmp(:)]));
            else
                stepAngle = skewAngle*(pi/180)/(Nz-1);
                p = zeros(mzptr.Nn*Nz,3);
                for i = 1:Nz
                    p((i-1)*mzptr.Nn+1:i*mzptr.Nn,1:3) = ...
                        ext_protate3z([mzptr.nodes,repmat(z(i),mzptr.Nn,1)],(i-1)*stepAngle);
                end
                ttmptr.addmz(mzName, TTMZPC(tmzpc_getExtrude(...
                    mzptr.cl,obj.elements(obj.ezi(:,mzptr.zi),1:3),...
                    mzptr.Nn,Nz-1),p));
            end
            ttmptr.mzs.(mzName).material = mzptr.material;
        end

        %% Index Operations
        % free boundary
        function y = getfbe(obj)
            y = find(obj.bedges);
        end
        function y = getfbn(obj)
            y = obj.getnIndexOnEdges(obj.getfbe);
        end
        % node index
        function y = getnIndexOnCircle(obj, c, r, tol)
            if nargin<4
                tol = obj.gleps;
            end
            y = [obj.nodes(:,1)-c(1),obj.nodes(:,2)-c(2)];
            y = sqrt(sum(y.^2,2));
            y = find(abs(y-r)<tol);
        end
        function y = getnIndexOnLine(obj, p1, p2, tol)
            if nargin<4
                tol = obj.gleps;
            end
            u = (p2-p1)/norm(p2-p1);
            pp1 = [obj.nodes(:,1)-p1(1),obj.nodes(:,2)-p1(2)];
            alpha = (pp1*u');
            y = sqrt(sum((pp1-alpha*u).^2,2));
            y = find(y<tol);
        end
        function y = getnIndexOnRay(obj,p1,p2)
            u = (p2-p1)/norm(p2-p1);
            pp1 = [obj.nodes(:,1)-p1(1),obj.nodes(:,2)-p1(2)];
            alpha = (pp1*u');
            y = sqrt(sum((pp1-alpha*u).^2,2));
            y = find(y<obj.gleps);
            index = alpha(y)<-obj.gleps;
            y = y(~index);
        end
        function y = getnIndexOnEdges(obj,eList)
            switch obj.etype
                case 'TL3'
                    y = obj.edges(eList,1:2);
                    y = unique(y(:));
                case 'TL6'
                    y = obj.edges(eList,1:2);
                    y = unique(y(:));
                    y = [y;eList+obj.Nn-size(obj.edges,1)];
            end
        end
        % edge index
        function y = geteIndexOnCircle(obj,varargin)
            y = obj.getnIndexOnCircle(varargin{:});
            y = bitand(ismember(obj.edges(:,1),y),ismember(obj.edges(:,2),y));
            y = find(y);
        end
        function y = geteIndexOnLine(obj,varargin)
            y = obj.getnIndexOnLine(varargin{:});
            y = bitand(ismember(obj.edges(:,1),y),ismember(obj.edges(:,2),y));
            y = find(y);
        end
        % periodic nodes
        function [km,ks] = splitRotate(obj,k,varargin)
            Nk = length(k);
            if rem(Nk,2)~=0
                error('number of input indices must be even.');
            end
            polarAngle = atan_02pi(obj.nodes(k,:));
            [~,index] = sort(polarAngle);
            k = k(index);
            pm = obj.nodes(k(1:Nk/2),:);
            ps = obj.nodes(k(Nk/2+1:Nk),:);
            km = k(1:Nk/2);
            ks = k(Nk/2+1:Nk);
            [Flag,index] = ismembertol(ps,pm,obj.gleps,'ByRow',true);
            if any(~Flag)
                error('These set of points does not form a set of peridic points.');
            end
            ks = ks(index);
        end
        function [km,ks] = splitPeriodic(obj,k,varargin)
            Nk = length(k);
            if rem(Nk,2)~=0
                error('number of input indices must be even.');
            end
            tp = obj.nodes(k,:);

            temp = 1;
            while ~isempty(tp)
                sp = tp(1,:);
                tp = tp(2:end,:);
                km(temp) = k(1);
                k = k(2:end);
                sp = ext_protate2(sp,varargin{1});
                index = find(sqrt(sum([tp(:,1)-sp(1),tp(:,2)-sp(2)].^2,2))<obj.gleps);
                if index
                    ks(temp) = k(index);
                    tp = tp([1:index-1,index+1:end],:);
                    k = k([1:index-1,index+1:end]);
                    temp = temp + 1;
                    continue
                end
                sp = ext_protate2(sp,-2*varargin{1});
                index = find(sqrt(sum([tp(:,1)-sp(1),tp(:,2)-sp(2)].^2,2))<obj.gleps);
                if index
                    ks(temp) = km(temp);
                    km(temp) = k(index);
                    tp = tp([1:index-1,index+1:end],:);
                    k = k([1:index-1,index+1:end]);
                    temp = temp + 1;
                    continue
                end
                error('These set of points does not form a set of peridic points.');
            end
        end
        %% Import and Export
        function read_g2d_txt(obj, Dir, MG)
            if nargin<3
                MG = 'MM';
            else
                MG = upper(strrep(MG,' ',''));
            end
            f = fopen(Dir,'r');
            clc
            tic
            while true
                str = fgetl(f);
                if isnumeric(str)
                    break
                end
                str = strsplit(str,' ');
                if strcmpi(str{1},'BeginFace:')
                    mzname = str{2};
                    str = strsplit(fgetl(f),' ');
                    Nloop = str2double(str{1});
                    f1 = cell(1,Nloop);
                    f2 = cell(1,Nloop);
                    v1 = cell(1,Nloop);
                    v2 = cell(1,Nloop);
                    for i = 1:Nloop
                        str = strsplit(fgetl(f),' ');
                        if strcmpi(str{1},'BeginLoop:')
                            % reading connectivity
                            str = fgetl(f);
                            Nf = str2double(str);
                            f1{i} = zeros(1,Nf);
                            f2{i} = zeros(1,Nf);
                            for j = 1:Nf
                                str = fgetl(f);
                                str = strsplit(str);
                                f1{i}(j) = str2double(str{1});
                                f2{i}(j) = str2double(str{2});
                            end
                            if i>1
                                f1{i} = f1{i} + length(v1{i-1});
                                f2{i} = f2{i} + length(v1{i-1});
                            end
                            % reading vertices
                            str = fgetl(f);
                            Nv = str2double(str);
                            v1{i} = zeros(1,Nv);
                            v2{i} = zeros(1,Nv);
                            for j = 1:Nv
                                str = fgetl(f);
                                str = strsplit(str);
                                v1{i}(j) = str2double(str{1});
                                v2{i}(j) = str2double(str{2});
                            end
                            str = fgetl(f);
                            str = strsplit(str);
                            if ~strcmp(str{1},'EndLoop.')
                                error('Wrong File Type.');
                            end
                        end
                    end
                    str = fgetl(f);
                    str = strsplit(str);
                    if ~strcmpi(str{1},'EndFace.')
                        error('Wrong File Type.');
                    end
                end
                % add to mesh database
                con = [cell2mat(f1);cell2mat(f2)]';
                ver = [cell2mat(v1);cell2mat(v2)]';
                [ver,~,ic] = uniquetol(ver,1e-6,'ByRows',true);
                con = ic(con);
                switch MG
                    case 'MM'
                        obj.addmz(mzname,MinimalMesh(con,ver));
                    case 'MG0'
                        obj.addmz(mzname,MeshGenerator0(con,ver));
                    case 'MG1'
                        obj.addmz(mzname,MeshGenerator1(con,ver));
                    otherwise
                        error('Mesh Generator does not exist ...');
                end
            end
            fclose(f);
            toc
        end
        function read_g2d_bin(obj, Dir, MG, varargin)
            if nargin<3
                MG = 'MM';
            end
            f = fopen(Dir,'r');
            clc
            tic
            while true
                str = fread(f,80,'*char');
                if isempty(str)
                    break;
                end
                str = strsplit(strtrim(str'));
                if strcmp(str{1},'BeginFace:')
                    mzname = str{2};
                    Nloop = fread(f,1,'uint32');
                    f1 = cell(1,Nloop);
                    f2 = cell(1,Nloop);
                    v1 = cell(1,Nloop);
                    v2 = cell(1,Nloop);
                    for i = 1:Nloop
                        str = fread(f,80,'*char');
                        str = strsplit(strtrim(str'));
                        if strcmp(str{1},'BeginLoop:')
                            % reading connectivity
                            Nf = fread(f,1,'uint32');
                            f1{i} = fread(f,Nf,'uint32')';
                            f2{i} = fread(f,Nf,'uint32')';
                            if i>1
                                f1{i} = f1{i} + length(v1{i-1});
                                f2{i} = f2{i} + length(v1{i-1});
                            end
                            % reading vertices
                            Nv = fread(f,1,'uint32');
                            v1{i} = fread(f,Nv,'double')';
                            v2{i} = fread(f,Nv,'double')';
                            str = fread(f,80,'*char');
                            str = strsplit(strtrim(str'));
                            if ~strcmp(str{1},'EndLoop')
                                error('Wrong File Type.');
                            end
                        end
                    end
                    str = fread(f,80,'*char');
                    str = strsplit(strtrim(str'));
                    if ~strcmp(str{1},'EndFace')
                        error('Wrong File Type.');
                    end
                end
                % add to mesh database
                con = [cell2mat(f1);cell2mat(f2)]';
                ver = [cell2mat(v1);cell2mat(v2)]';
                [ver,~,ic] = uniquetol(ver,1e-6,'ByRows',true);
                con = ic(con);
                switch MG
                    case 'MM'
                        obj.addmz(mzname,MinimalMesh(con,ver, varargin{:}));
                    case 'MG0'
                        obj.addmz(mzname,MeshGenerator0(con,ver, varargin{:}));
                    case 'MG1'
                        obj.addmz(mzname,MeshGenerator1(con,ver, varargin{:}));
                    otherwise
                        error('Invalid mesh generator name.');
                end
            end
            fclose(f);
            toc
        end
        function read_stlf_bin(obj, Dir)
            if ~isdir(Dir)
                error('Directory was not found.');
            end
            stlFile = ls(Dir);
            stlFile = stlFile(3:end,:);
            for i = 1:size(stlFile,1)
                obj.read_stl_bin([Dir,'\',stlFile(i,:)]);
            end
        end
        function read_stl_bin(obj, Dir)
            % openning STL file
            f = fopen(Dir,'r');
            if f<0
                error('Can not open file.');
            end
            % zone name
            mzname = fread(f,80,'*char');
            mzname = strrep(mzname', ' ', '');
            % number of elements
            Nt = fread(f,1,'uint32');
            % trianle points
            p = zeros(9,Nt);
            % loop over elements
            for i = 1:Nt
                % reading normal vector
                fread(f,3,'single');
                % reading vertices
                p(:,i) = fread(f,9,'single');
                % reading color
                fread(f,1,'uint16');
            end
            % closing STL file
            fclose(f);
            % generation of new mesh zone
            % nodes
            p = reshape(p,3,[]);
            p = p';
            % elements
            e = 1:3*Nt;
            e = reshape(e,3,[]);
            e = e';
            % unification
            [p,~,index] = uniquetol(p,1e-6,'ByRows',true);
            % reindex
            e = index(e);
            % adding mesh zone
            obj.addmz(mzname,TMZPC(e,p(:,1:2)));
        end
        function write_stl_bin(obj)
            % openning STL file
            if exist('STL_Files','file')
                rmdir('STL_Files','s')
            end
            mkdir('STL_Files')
            mznames = fieldnames(obj.mzs);
            for i = 1:numel(mznames)
                % zone name
                mzname = mznames{i};
                f = fopen(['STL_Files\', mzname, '.STL'],'w');
                % writing zone name
                fwrite(f,[mzname,repmat(' ',1,80-length(mzname))],'*char');
                % number of elements
                fwrite(f,obj.mzs.(mzname).Ne,'uint32');
                % loop over elements
                for j = 1:obj.mzs.(mzname).Ne
                    % writing normal vector
                    fwrite(f,[0,0,1],'single');
                    % writing vertices
                    ti = obj.mzs.(mzname).cl(j,:);
                    fwrite(f,[obj.mzs.(mzname).nodes(ti(1),:),0]','single');
                    fwrite(f,[obj.mzs.(mzname).nodes(ti(2),:),0]','single');
                    fwrite(f,[obj.mzs.(mzname).nodes(ti(3),:),0]','single');
                    % writing color
                    fwrite(f,1,'uint16');
                end
                % closing STL file
                fclose(f);
            end
        end
        function write_m2d_bin(obj, Dir)
            % openning STL foler
            if nargin<2
                if exist('m2d_Files','file')
                    rmdir('m2d_Files','s');
                end
                mkdir('m2d_Files');
                Dir = 'm2d_Files';
            end
            % loop over mesh zones
            mznames = fieldnames(obj.mzs);
            for i = 1:numel(mznames)
                % zone name
                mzname = mznames{i};
                f = fopen([Dir,'\', mzname, '.m2d'],'w');
                % writing number of nodes and elements
                fwrite(f,uint32(obj.mzs.(mzname).Nn),'uint32');
                fwrite(f,uint32(obj.mzs.(mzname).Ne),'uint32');
                % writing nodes
                fwrite(f,obj.mzs.(mzname).nodes(:),'double');
                % writing elements
                fwrite(f,obj.mzs.(mzname).cl(:),'double');
                % closing m2d file
                fclose(f);
            end
        end
        function read_m2d_bin(obj, Dir)
            % openning m2d file
            f = fopen(Dir,'r');
            Nnodes = fread(f,1,'uint32');
            Nelements = fread(f,1,'uint32');
            p = fread(f,2*Nnodes,'double');
            e = fread(f,3*Nelements,'double');
            fclose(f);
            p = reshape(p,[],2);
            e = reshape(e,[],3);
            % adding mesh zone
            mzname = strsplit(Dir, '\');
            mzname = strsplit(mzname{end}, '.');
            obj.addmz(mzname{1},TMZPC(e,p));
        end
        function read_m2df_bin(obj, Dir)
            if ~isdir(Dir)
                error('Directory was not found.');
            end
            stlFile = ls(Dir);
            stlFile = stlFile(3:end,:);
            for i = 1:size(stlFile,1)
                obj.read_m2d_bin([Dir,'\',stlFile(i,:)]);
            end
        end
        function read_msh_bin(obj, Dir)
            % openning msh file
            f = fopen(Dir,'r');
            if f<0
                error('Can not open file.');
            end
            % allocating space for mesh zone elements
            e = cell(1,20);
            for i = 1:numel(e)
                e{i} = zeros(3,[]);
            end
            % loop for read
            while true
                str =  fgetl(f);
                if strcmp(str,'$Nodes')
                    Nnodes = fgetl(f);
                    Nnodes = str2double(Nnodes);
                    p = zeros(3,Nnodes);
                    for i = 1:Nnodes
                        fread(f,1,'int32');
                        p(:,i) = fread(f,3,'double');
                    end
                end
                if strcmp(str,'$Elements')
                    Nnodes = fgetl(f);
                    Nnodes = str2double(Nnodes);
                    for i = 1:Nnodes
                        eIndex = fread(f,1,'int32');
                        if eIndex == 15
                            Nelements = fread(f,1,'int32');
                            Nt = fread(f,1,'int32');
                            for j = 1:Nelements
                                fread(f,1,'int32');
                                fread(f,Nt,'int32');
                                fread(f,1,'int32');
                            end
                        elseif eIndex == 1
                            Nelements = fread(f,1,'int32');
                            Nt = fread(f,1,'int32');
                            for j = 1:Nelements
                                fread(f,1,'int32');
                                fread(f,Nt,'int32');
                                fread(f,2,'int32');
                            end
                        elseif eIndex == 2
                            Nelements = fread(f,1,'int32');
                            Nt = fread(f,1,'int32');
                            for j = 1:Nelements
                                fread(f,1,'int32');
                                fread(f,Nt-1,'int32');
                                zIndex = fread(f,1,'int32');
                                e{zIndex}(:,end+1) = fread(f,3,'int32');
                            end
                        end
                    end
                end
                if strcmp(str,'$EndElements'), break; end
            end
            fclose(f);
            % adding mesh zone
            e = e(1:zIndex);
            p = p(1:2,:)';
            for i = 1:zIndex
                xtmp = p(:,1);
                ytmp = p(:,2);
                xtmp = xtmp(e{i});
                ytmp = ytmp(e{i});
                [ptmp, ~, index] = uniquetol([xtmp(:),ytmp(:)], 1e-6, 'ByRows', true);
                etmp = 1:3*size(e{i},2);
                etmp = etmp(index);
                etmp = reshape(etmp,3,[])';
                obj.addmz(TMZPC(etmp,ptmp));
            end
        end
        function aux_unify(obj, mzNameCommonPart)
            mzNames = obj.getMeshZoneNames;
            idx = ~cellfun('isempty', regexp(mzNames, "^"+mzNameCommonPart+"\d+$"));
            if isempty(idx)
                error('No mesh zone was detected for this common part string.');
            else
                obj.joinMeshZones(mzNameCommonPart, mzNames(idx));
            end
        end

        % operation sequence: copy mirror -> copy rotate
        % there are two outputs
        % output is a list containing the name of new mesh zones
        function varargout = aux_cmcr(obj, mzName, mirrorAxis, Ncopy, rotateAngle, xc, yc)

            if nargin < 5
                rotateAngle = 2*pi/Ncopy;
            end

            if nargin < 6
                xc = 0;
                yc = 0;
            end

            % change the name of current mesh zone
            obj.changeMeshZoneName(mzName, [mzName, '11']);

            % mirror mesh zone with respect to x axis
            obj.cmmz([mzName, '21'], [mzName, '11'], mirrorAxis);

            for i = 2:Ncopy
                obj.crmz([mzName, '1', num2str(i)], [mzName, '1', num2str(i-1)], rotateAngle, [xc, yc]);
                obj.crmz([mzName, '2', num2str(i)], [mzName, '2', num2str(i-1)], rotateAngle, [xc, yc]);
            end

            if nargout == 1

                newMeshZoneNames = strings(1,2*Ncopy);
                for i = 1:Ncopy
                    newMeshZoneNames(2*i-1) = mzName + "1" + num2str(i);
                    newMeshZoneNames(2*i) = mzName + "2" + num2str(i);
                end
                varargout{1} = newMeshZoneNames;

            elseif nargout > 1
                error('Too many output arguments');
            end

            % change states
            obj.makeFalse_isGlobalMeshGenerated;

        end

        function aux_cmjcrj(obj, mzName, varargin)
            newNames = obj.aux_cmcr(mzName, varargin{:});
            obj.joinMeshZones(mzName, newNames);
        end

        function aux_cmxjcrj(obj, mzName, varargin)
            mzName = char(mzName);
            newNames = obj.aux_cmcr(mzName, [1,0], varargin{:});
            obj.joinMeshZones(mzName, newNames);
        end

        function aux_cmxjcr(obj, mzName, varargin)
            mzName = char(mzName);
            newNames = obj.aux_cmcr(mzName, [1,0], varargin{:});
            for i = 1:(numel(newNames)/2)
                newName = char(mzName + string(i));
                obj.joinMeshZones(newName, newNames(2*i-1:2*i));
            end
        end
        
    end
    methods (Access = private)
        function makeFalse_isGlobalMeshGenerated(obj)
            obj.isGlobalMeshGenerated = false;
            obj.isd2ElementsGenerated = false;
            obj.isd3ElementsGenerated = false;
            obj.isJITEvaluated = false;
            obj.isKeFe_TL3_Evaluated = false;
            obj.isKeFe_TL6_Evaluated = false;
            obj.isKexy4Fe_TL3_Evaluated = false;
            obj.isKeFe_TL6_cte_Evaluated = false;
            obj.isKe_TL3_Evaluated = false;
            obj.isFe_TL3_Evaluated = false;
            obj.isMe_TL3_Evaluated = false;
            obj.isKexy1_TL3_Evaluated = false;
        end
    end
    
end

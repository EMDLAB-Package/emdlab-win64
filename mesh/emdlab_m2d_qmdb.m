% EMDLAB: Electrical Machines Design Laboratory
% 2D quadrilateral mesh data base class

classdef emdlab_m2d_qmdb < handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_mdb_cp & emdlab_m2d_xmdb

    properties (Dependent = true)

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
            % set default element type
            obj.etype = 'QL4';
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
            %             obj.mzs.(mzName).color = rand(1,3);

            % changing states
            obj.clearGlobalMeshGenerationFlag;

        end

        function delete(obj)

            mzNames = fieldnames(obj.mzs);

            for i = 1:numel(mzNames)
                delete(obj.mzs.(mzNames{i}));
            end

        end

        %% FEM Preparation
        function ggmesh(obj)
            % generate global mesh

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
            obj.gea = zeros(1,obj.Ne);
            emdlab_m2d_qmdbc_evalgea(obj.cl, obj.nodes, obj.gea);
            obj.gea = abs(obj.gea);

            % change states
            obj.isGlobalMeshGenerated = true;

            % settig element type
            obj.etype = 'QL4';

        end

        function evalezi(obj)
            % evaluate element material & zone index

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
                obj.clearGlobalMeshGenerationFlag;
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
                obj.clearGlobalMeshGenerationFlag;
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

            % initialize nbs matrix
            obj.nbs = zeros(obj.Ne,4);

            % edge length
            obj.edgeLength = sqrt(sum((obj.nodes(obj.edges(:,1),:) - obj.nodes(obj.edges(:,2),:)).^2, 2));
            obj.el = obj.edgeLength(abs(obj.elements(:,1:4)));

            % edge unit vectors
            obj.uEdges = obj.nodes(obj.edges(:,2),:) - obj.nodes(obj.edges(:,1),:);
            obj.uEdges = obj.uEdges ./ vecnorm(obj.uEdges, 2, 2);
            obj.nEdges = [obj.uEdges(:,2), -obj.uEdges(:,1)];

            % edge element
            obj.edges = [obj.edges, zeros(size(obj.edges,1), 6)];
            emdlab_m2d_qmdbc_evalee(obj.edges, obj.elements, obj.nbs);

        end

        function strefine(obj)
            mznames = fieldnames(obj.mzs);
            for i = 1:numel(mznames)
                obj.mzs.(mznames{i}).strefine;
            end
            obj.clearGlobalMeshGenerationFlag;
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

        %% Tools: some operations on mesh zones
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

        function y = getCenterOfElements(obj)
            obj.ggmesh;
            y = (obj.nodes(obj.cl(:, 1), :) + obj.nodes(obj.cl(:, 2), :) + obj.nodes(obj.cl(:, 3), :) + obj.nodes(obj.cl(:, 4), :)) / 4;
        end

        function y = getCenterOfEdges(obj)
            obj.ggmesh;
            y = (obj.nodes(obj.edges(:, 1), :) + obj.nodes(obj.edges(:, 2), :)) / 2;
        end

        %% Index Finding Functions
        function y = getfbn(obj)
            y = obj.getnIndexOnEdges(obj.getfbe);
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

        %% Named Selections
        % edges
        function name = checkEdgeNamedSelectionExistence(obj, name)

            name = erase(name, ' ');
            if ~isfield(obj.edgeNamedSelections, name)
                error('Specified edge named selection does not exist.');
            end

        end

        function name = checkEdgeNamedSelectionNonExistence(obj, name)

            name = erase(name, ' ');
            if isfield(obj.edgeNamedSelections, name)
                error('Specified edge named selection already exist.');
            end

        end

        function addEdgeNamedSelection(obj, name, indices)

            name = erase(name, ' ');
            name = obj.checkEdgeNamedSelectionNonExistence(name);
            if isempty(indices)
                error('Indices matreix cannot be empty.');
            end

            obj.edgeNamedSelections.(name) = indices;

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

        %% Auxiliary Functions
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
            obj.clearGlobalMeshGenerationFlag;

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

        %% Solver Related Functions (FVM)
        function [indexI,indexJ,value,sourceVector] = buildGMatrixSVector(obj, resistances)

            obj.ggmesh;
            indexI = zeros(obj.Ne,5);
            indexJ = zeros(obj.Ne,5);
            value = zeros(obj.Ne,5);
            sourceVector = zeros(obj.Ne,1);

            % loop over elements to construct G_matrix
            for i = 1:obj.Ne

                % constructing i-th row
                indexI(i,:) = i;

                % conductance between cell P and 1 nb
                if obj.nbs(i,1)
                    if ismember(abs(obj.elements(i,1)),abs(obj.elements(obj.nbs(i,1),[1,3])))
                        Rp1 = resistances(i,1)/2 + resistances(obj.nbs(i,1),1)/2;
                    else
                        Rp1 = resistances(i,1)/2 + resistances(obj.nbs(i,1),2)/2;
                    end
                    indexJ(i,2) = obj.nbs(i,1);
                    value(i,2) = -1/Rp1;
                else
                    Rp1 = inf;
                    indexJ(i,2) = i;
                end

                % conductance between cell P and 2 nb
                if obj.nbs(i,2)
                    if ismember(abs(obj.elements(i,2)),abs(obj.elements(obj.nbs(i,2),[1,3])))
                        Rp2 = resistances(i,2)/2 + resistances(obj.nbs(i,2),1)/2;
                    else
                        Rp2 = resistances(i,2)/2 + resistances(obj.nbs(i,2),2)/2;
                    end
                    indexJ(i,3) = obj.nbs(i,2);
                    value(i,3) = -1/Rp2;
                else
                    Rp2 = inf;
                    indexJ(i,3) = i;
                end

                % conductance between cell P and 3 nb
                if obj.nbs(i,3)
                    if ismember(abs(obj.elements(i,3)),abs(obj.elements(obj.nbs(i,3),[1,3])))
                        Rp3 = resistances(i,1)/2 + resistances(obj.nbs(i,3),1)/2;
                    else
                        Rp3 = resistances(i,1)/2 + resistances(obj.nbs(i,3),2)/2;
                    end
                    indexJ(i,4) = obj.nbs(i,3);
                    value(i,4) = -1/Rp3;
                else
                    Rp3 = inf;
                    indexJ(i,4) = i;
                end

                % conductance between cell P and 4 nb
                if obj.nbs(i,4)
                    if ismember(abs(obj.elements(i,4)),abs(obj.elements(obj.nbs(i,4),[1,3])))
                        Rp4 = resistances(i,2)/2 + resistances(obj.nbs(i,4),1)/2;
                    else
                        Rp4 = resistances(i,2)/2 + resistances(obj.nbs(i,4),2)/2;
                    end
                    indexJ(i,5) = obj.nbs(i,4);
                    value(i,5) = -1/Rp4;
                else
                    Rp4 = inf;
                    indexJ(i,5) = i;
                end

                indexJ(i,1) = i;
                value(i,1) = 1/Rp1 + 1/Rp2 + 1/Rp3 + 1/Rp4;

            end
        end

        function [value, sourceVector] = applyFixedTemperatureBC(obj, idx, value, sourceVector, resistances, TValue)

            if iscolumn(idx), idx = idx'; end
            edgeCenters = obj.getCenterOfEdges;

            for index = idx

                if isa(TValue,'function_handle')
                    tvbc = TValue(edgeCenters(index,1),edgeCenters(index,2));
                else
                    tvbc = TValue;
                end

                if obj.edges(index,5)

                    eIndex = obj.edges(index,5);
                    if ismember(obj.edges(index,6),[1,3])
                        value(eIndex,1) = value(eIndex,1) + 2/resistances(eIndex,1);
                        sourceVector(eIndex) = sourceVector(eIndex) + tvbc*2/resistances(eIndex,1);
                    else
                        value(eIndex,1) = value(eIndex,1) + 2/resistances(eIndex,2);
                        sourceVector(eIndex) = sourceVector(eIndex) + tvbc*2/resistances(eIndex,2);
                    end

                else

                    eIndex = obj.edges(index,7);
                    if ismember(obj.edges(index,8),[1,3])
                        value(eIndex,1) = value(eIndex,1) + 2/resistances(eIndex,1);
                        sourceVector(eIndex) = sourceVector(eIndex) + tvbc*2/resistances(eIndex,1);
                    else
                        value(eIndex,1) = value(eIndex,1) + 2/resistances(eIndex,2);
                        sourceVector(eIndex) = sourceVector(eIndex) + tvbc*2/resistances(eIndex,2);
                    end

                end

            end

        end

        function [value, sourceVector] = applyConvectionBC(obj, idx, value, sourceVector,resistances,z,Tinf, HValue,unit)

            if iscolumn(idx), idx = idx'; end
            edgeCenters = obj.getCenterOfEdges;

            for index = idx

                if isa(HValue,'function_handle')
                    hvbc = HValue(edgeCenters(index,1),edgeCenters(index,2));
                else
                    hvbc = HValue;
                end

                if obj.edges(index,5)

                    eIndex = obj.edges(index,5);
                    if ismember(obj.edges(index,6),[1,3])
                        Gconv = 1/(resistances(eIndex,1)/2 + 1/(hvbc*z*unit*obj.el(eIndex,obj.edges(index,6))));
                    else
                        Gconv = 1/(resistances(eIndex,2)/2 + 1/(hvbc*z*unit*obj.el(eIndex,obj.edges(index,6))));
                    end

                else

                    eIndex = obj.edges(index,7);
                    if ismember(obj.edges(index,8),[1,3])
                        Gconv = 1/(resistances(eIndex,1)/2 + 1/(hvbc*z*unit*obj.el(eIndex,obj.edges(index,8))));
                    else
                        Gconv = 1/(resistances(eIndex,2)/2 + 1/(hvbc*z*unit*obj.el(eIndex,obj.edges(index,8))));
                    end

                end

                value(eIndex,1) = value(eIndex,1) + Gconv;
                sourceVector(eIndex) = sourceVector(eIndex)  + Tinf*Gconv;

            end

        end

        function [value, sourceVector] = applyHeatFluxBC(obj, idx, value, sourceVector, HFValue,unit)

            if iscolumn(idx), idx = idx'; end
            edgeCenters = obj.getCenterOfEdges;

            for index = idx

                if isa(HFValue,'function_handle')
                    hfbc = HFValue(edgeCenters(index,1),edgeCenters(index,2));
                else
                    hfbc = HFValue;
                end

                if obj.edges(index,5)

                    eIndex = obj.edges(index,5);
                    sourceVector(eIndex) = sourceVector(eIndex) + hfbc*obj.el(eIndex,obj.edges(index,6))*unit*0.01;

                else

                    eIndex = obj.edges(index,7);
                    sourceVector(eIndex) = sourceVector(eIndex) + hfbc*obj.el(eIndex,obj.edges(index,8))*unit*0.01;

                end

            end

        end

        function sourceVector = applyIHG(obj, mzName, sourceVector,z, IHGValue, unit)

            idx = 1:obj.Ne;
            idx = idx(obj.ezi(:,obj.mzs.(mzName).zi));

            for index = idx

                if isa(IHGValue,'function_handle')
                    ihgbc = IHGValue(edgeCenters(index,1),edgeCenters(index,2));
                else
                    ihgbc = IHGValue;
                end

                sourceVector(index) = sourceVector(index) + ihgbc;

            end

        end

        function varargout = showAverageTemperature(obj, T, N, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;
            if nargin<3, N=10; end

            patch('faces', obj.cl, 'Vertices', obj.nodes, 'FaceVertexCData',T,'FaceColor','flat', 'edgecolor', 'none');
            colormap(jet(N));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Average Temperature';

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

    end

    methods (Access = private)

        function clearGlobalMeshGenerationFlag(obj)

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

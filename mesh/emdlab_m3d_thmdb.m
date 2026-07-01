% EMDLAB: Electrical Machines Design Laboratory
% tetrahedral mesh data base

classdef emdlab_m3d_thmdb < handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_mdb_cp

    properties (SetAccess = private)

        % mesh nodes: [x,y,z]
        nodes (:,3) double;

        % mesh connectivity list
        cl (:,:) double;

        % mesh elements: [facet1, facet2, facet3, facet4, zone index]
        elements (:,5) double;

        % unique edges (:,8): [index of node1, index of node2, ]
        edges

        % Unique Mesh Facets
        facets

        % Boundary Facets
        bfacets

        % jacobian inverse transpose
        JIT (9,:) double;

        % element zone index
        ezi (:,:) logical;

        % global elements volume
        gev (1,:) double;

        % auxiliary stored matricies
        mtcs (1,1) struct;

        % named selections
        facetNamedSelections (1,1) struct;

        % flag to print the elapsed times
        printFlag (1,1) logical = true;

        % element type
        etype (1,:) char = 'TTL4';

    end

    properties (Dependent = true)

        % Number of nodes
        Nn (1,1) double;

        % Number of elements
        Ne (1,1) double;

        % flags for element type
        isTTL4 (1,1) logical;
        isTTL10 (1,1) logical;

    end

    properties (Access = private)

        % states
        isd2ElementsGenerated (1,1) logical = false;
        isd3ElementsGenerated (1,1) logical = false;
        isJITEvaluated = false;
        isKeMeFe_TTL4_Evaluated (1,1) logical = false;

    end

    methods
        %% constructor and destructor
        function obj = emdlab_m3d_thmdb(varargin)
            % add default material
            obj.addMaterial('air');
        end

        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end

        function y = get.Ne(obj)
            y = size(obj.cl, 1);
        end

        function setPrintFlag(obj, newValue)
            obj.printFlag = newValue;
        end

        function delete(obj)

            mzNames = fieldnames(obj.mzs);

            for i = 1:numel(mzNames)
                delete(obj.mzs.(mzNames{i}));
            end

        end

        function y = getMeshZoneNames(obj)
            y = string(fieldnames(obj.mzs))';
        end

        function y = get.isTTL4(obj)
            y = strcmpi(obj.etype, 'TTL4');
        end

        function y = get.isTTL10(obj)
            y = strcmpi(obj.etype, 'TTL10');
        end

        %% FEM preparation
        % generate global mesh
        function ggmesh(obj, mzFlag)

            if nargin<2, mzFlag = false; end

            % check states
            if obj.isGlobalMeshGenerated, return; end

            % generation of initial mesh
            Nn_tmp = 0;
            Ne_tmp = 0;
            mzNames = fieldnames(obj.mzs);

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                Nn_tmp = Nn_tmp + mzptr.Nn;
                Ne_tmp = Ne_tmp + mzptr.Ne;
            end

            % initialization of nodes and elements
            obj.nodes = zeros(Nn_tmp, 3);
            obj.cl = zeros(Ne_tmp, 4);
            obj.elements = zeros(Ne_tmp, 5);
            nindex = 0;
            eindex = 0;

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                % insertion of nodes
                obj.nodes(1 + nindex:mzptr.Nn + nindex, :) = mzptr.nodes;
                % insertion of elements
                obj.cl(1 + eindex:mzptr.Ne + eindex, :) = mzptr.cl + nindex;
                % specefying zone index
                obj.elements(1 + eindex:mzptr.Ne + eindex, 5) = i;
                mzptr.zi = i;
                nindex = nindex + mzptr.Nn;
                eindex = eindex + mzptr.Ne;
            end

            [obj.nodes, ~, ic] = uniquetol(obj.nodes, obj.gleps, 'ByRows', true);
            obj.cl = ic(obj.cl);
            % setting l2g
            nindex = 0;

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                mzptr.l2g = ic(nindex + 1:nindex + mzptr.Nn);
                nindex = nindex + mzptr.Nn;
            end

            obj.setdata;
            obj.evalezi;

            if mzFlag

                for i = 1:obj.Nmzs
                    mzptr =  obj.mzs.(mzNames{i});
                    mzptr.props.cl = obj.cl(obj.ezi(:,mzptr.zi),:);
                end

            end

            % change states
            obj.isGlobalMeshGenerated = true;

            % settig element type
            obj.etype = 'TTL4';

        end

        % evaluate element zone index
        function evalezi(obj)

            obj.ezi = false(obj.Ne, obj.Nmzs);

            for i = 1:obj.Nmzs
                obj.ezi(:, i) = obj.elements(:, 5) == i;
            end

        end

        % evaluate jacobian inverse transpose matrix
        function evalJIT(obj, mzFlag)

            if nargin<2, mzFlag = false; end

            % check states
            if obj.isJITEvaluated, return; end
            if obj.printFlag
                tic, disp('-------------------------------------------------------');
            end

            % prerequisite
            obj.ggmesh(mzFlag);

            % evaluation of jacobian inverse transpose using d1 element data
            % JIT is a [9 x Ne] matrix [J11;J21;J31;J12;J22;J32;J13;J23;J33]
            [obj.JIT, obj.gev] = emdlab_m3d_ttl4_evalJIT(obj.cl(:,1:4), obj.nodes);

            if mzFlag

                mzNames = fieldnames(obj.mzs);
                for i = 1:obj.Nmzs
                    mzptr = obj.mzs.(mzNames{i});
                    mzptr.props.JIT = obj.JIT(:,obj.ezi(:,mzptr.zi));
                    mzptr.props.gev = obj.gev(obj.ezi(:,mzptr.zi));
                end

            end

            % change states
            obj.isJITEvaluated = true;
            if obj.printFlag
                disp('Evaluation of JIT completed.');
                toc, disp('-------------------------------------------------------');
            end

        end

        % tetrahedral elements: first order
        function emdlab_m3d_ttl4_evalKeMeFeNodal(obj)

            % prerequests
            obj.ggmesh;
            obj.evalJIT;

            tic, disp('-------------------------------------------------------');
            if obj.isKeMeFe_TTL4_Evaluated, return; end

            % calculation of stiffness matrix, mass matrix and force vector
            [obj.mtcs.Kex, obj.mtcs.Key, obj.mtcs.Kez, obj.mtcs.Me, obj.mtcs.Fe] = emdlab_m3d_ttl4_evalKeMeFeNodal(obj.cl, obj.nodes);

            disp('Calculation of [Ke], [Me], and [Fe] completed.');
            toc, disp('-------------------------------------------------------');

            % change states
            obj.isKeMeFe_TTL4_Evaluated = true;

        end

        function evalKeFe_TTL4(obj)
            if obj.isKeFe_TTL4_Evaluated, return; end
            obj.ggmesh;
            obj.evalJIT;
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic, disp('-------------------------------------------------------');
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke = zeros(10, xNe);
            temp = 0;

            for i = 1:4

                for j = 1:i
                    temp = temp + 1;
                    obj.Ke(temp, :) = ...
                        obj.gphix(i, :) .* obj.gphix(j, :) + ...
                        obj.gphiy(i, :) .* obj.gphiy(j, :) + ...
                        obj.gphiz(i, :) .* obj.gphiz(j, :);
                end

            end

            % multiplying by triangle areas
            obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gev);
            % elemental force matrix
            obj.Fe = repmat((obj.gev / 4), 4, 1);
            disp('Evaluation of [Ke] and [Fe] completed.');
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKeFe_TTL4_Evaluated = true;
        end

        function evalKeFe(obj, etype)
            xNe = obj.Ne;
            obj.etype = etype;
            % getting data of reference element
            tic
            % evaluation of jacobian inverse using d1 element data
            edata = getedata('TTL4');
            % x, y and z coordinate of points
            xp = obj.nodes(:, 1);
            yp = obj.nodes(:, 2);
            zp = obj.nodes(:, 3);
            % point coordinate of each triangle nodes
            xp = xp(obj.cl');
            yp = yp(obj.cl');
            zp = zp(obj.cl');
            % evaluation of a, b and c coefficient
            acoefs = edata.M \ xp;
            bcoefs = edata.M \ yp;
            ccoefs = edata.M \ zp;
            % setting elements volume
            obj.gev = zeros(1, obj.Ne);
            mznames = fieldnames(obj.mzs);

            for i = 1:obj.Nmzs
                obj.gev(obj.ezi(:, obj.mzs.(mznames{i}).zi)) = ...
                    obj.mzs.(mznames{i}).ev;
            end

            % Jacobiab transpose matrix of each elements
            obj.JT = [acoefs(2, :); acoefs(3, :); acoefs(4, :); ...
                bcoefs(2, :); bcoefs(3, :); bcoefs(4, :); ...
                ccoefs(2, :); ccoefs(3, :); ccoefs(4, :); ];
            [i, j] = getij(3, xNe);
            obj.JT = sparse(i, j, obj.JT(:));
            % getting specefied edata
            edata = getedata(obj.etype);

            switch obj.etype
                case 'TTL4'
                    gphi = obj.JT \ repmat(edata.GG, xNe, 1);
                    % inserting in gradx and grady phi i
                    obj.gphix = (gphi(1:3:end, :))';
                    obj.gphiy = (gphi(2:3:end, :))';
                    obj.gphiz = (gphi(3:3:end, :))';
                    % evaluation of Ke
                    obj.Ke = zeros(10, xNe);
                    temp = 0;

                    for i = 1:4

                        for j = 1:i
                            temp = temp + 1;
                            obj.Ke(temp, :) = ...
                                obj.gphix(i, :) .* obj.gphix(j, :) + ...
                                obj.gphiy(i, :) .* obj.gphiy(j, :) + ...
                                obj.gphiz(i, :) .* obj.gphiz(j, :);
                        end

                    end

                    % multiplying by triangle areas
                    obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gev);
                    % elemental force matrix
                    obj.Fe = repmat((obj.gev / 4), 4, 1);
                otherwise
                    error('Element type does not defined.');
            end

            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
        end

        function evalKeFeTTL4_Nodal(obj)
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke.K11 = zeros(10, xNe);
            obj.Ke.K22 = zeros(10, xNe);
            obj.Ke.K33 = zeros(10, xNe);
            obj.Ke.K21 = zeros(16, xNe);
            obj.Ke.K31 = zeros(16, xNe);
            obj.Ke.K32 = zeros(16, xNe);
            temp = 0;

            for i = 1:4

                for j = 1:i
                    temp = temp + 1;
                    obj.Ke.K11(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :) + obj.gphiz(i, :) .* obj.gphiz(j, :);
                    obj.Ke.K22(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :) + obj.gphiz(i, :) .* obj.gphiz(j, :);
                    obj.Ke.K33(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :) + obj.gphix(i, :) .* obj.gphix(j, :);
                end

            end

            temp = 0;

            for i = 1:4

                for j = 1:4
                    temp = temp + 1;
                    obj.Ke.K21(temp, :) = -1 * (obj.gphix(i, :) .* obj.gphiy(j, :));
                    obj.Ke.K31(temp, :) = -1 * (obj.gphix(i, :) .* obj.gphiz(j, :));
                    obj.Ke.K32(temp, :) = -1 * (obj.gphiy(i, :) .* obj.gphiz(j, :));
                end

            end

            % multiplying by triangle areas
            obj.Ke.K11 = obj.Ke.K11 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K22 = obj.Ke.K22 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K33 = obj.Ke.K33 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K21 = obj.Ke.K21 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K31 = obj.Ke.K31 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K32 = obj.Ke.K32 * sparse(1:xNe, 1:xNe, obj.gev);
            % elemental force matrix
            obj.Fe = repmat((obj.gev / 4), 4, 1);
            disp('Evaluation of [Ke] and [Fe] completed.')
            toc, disp('-------------------------------------------------------');
        end

        function evalKexyz1_TTL4(obj)
            if obj.isKexyz1_TTL4_Evaluated, return; end
            obj.ggmesh;
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke = zeros(10, obj.Ne);
            temp = 0;

            for i = 1:4

                for j = 1:i
                    temp = temp + 1;
                    obj.Ke(temp, :) = (obj.gphix(i, :) .* obj.gphix(j, :) + ...
                        obj.gphiy(i, :) .* obj.gphiy(j, :) + ...
                        obj.gphiz(i, :) .* obj.gphiz(j, :)) .* obj.gev;
                end

            end

            disp('Evaluation of [Ke] and [Fe] completed.')
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexyz1_TTL4_Evaluated = true;
        end

        function evalKexyz3_TTL4(obj)
            if obj.isKexyz3_TTL4_Evaluated, return; end
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke.xx = zeros(10, obj.Ne);
            obj.Ke.yy = zeros(10, obj.Ne);
            obj.Ke.zz = zeros(10, obj.Ne);
            temp = 0;

            for i = 1:4

                for j = 1:i
                    temp = temp + 1;
                    obj.Ke.xx(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :) .* obj.gev;
                    obj.Ke.yy(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :) .* obj.gev;
                    obj.Ke.zz(temp, :) = obj.gphiz(i, :) .* obj.gphiz(j, :) .* obj.gev;
                end

            end

            disp('Evaluation of [Ke] and [Fe] completed.')
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexyz3_TTL4_Evaluated = true;
        end

        function evalKexyz9FeTTL4(obj)
            if obj.isKexyz9Fe_Evaluated, return; end
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke.Kxx = zeros(10, xNe);
            obj.Ke.Kyy = zeros(10, xNe);
            obj.Ke.Kzz = zeros(10, xNe);
            obj.Ke.Kyx = zeros(16, xNe);
            obj.Ke.Kzx = zeros(16, xNe);
            obj.Ke.Kzy = zeros(16, xNe);
            temp = 0;

            for i = 1:4

                for j = 1:i
                    temp = temp + 1;
                    obj.Ke.Kxx(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :) .* obj.gev;
                    obj.Ke.Kyy(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :) .* obj.gev;
                    obj.Ke.Kzz(temp, :) = obj.gphiz(i, :) .* obj.gphiz(j, :) .* obj.gev;
                end

            end

            temp = 0;

            for i = 1:4

                for j = 1:4
                    temp = temp + 1;
                    obj.Ke.Kyx(temp, :) = obj.gphix(i, :) .* obj.gphiy(j, :) .* obj.gev;
                    obj.Ke.Kzx(temp, :) = obj.gphix(i, :) .* obj.gphiz(j, :) .* obj.gev;
                    obj.Ke.Kzy(temp, :) = obj.gphiy(i, :) .* obj.gphiz(j, :) .* obj.gev;
                end

            end

            % elemental force matrix
            obj.Fe = repmat((obj.gev / 4), 4, 1);
            disp('Evaluation of [Ke] and [Fe] completed.')
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexyz9Fe_Evaluated = true;
        end

        function obj = evalKe6Fe(obj, etype)
            xNe = obj.Ne;
            obj.etype = etype;
            % getting data of reference element
            tic
            % evaluation of jacobian inverse using d1 element data
            edata = getedata('TTL4');
            % x, y and z coordinate of points
            xp = obj.nodes(:, 1);
            yp = obj.nodes(:, 2);
            zp = obj.nodes(:, 3);
            % point coordinate of each triangle nodes
            xp = xp(obj.cl');
            yp = yp(obj.cl');
            zp = zp(obj.cl');
            % evaluation of a, b and c coefficient
            acoefs = edata.M \ xp;
            bcoefs = edata.M \ yp;
            ccoefs = edata.M \ zp;
            % setting elements volume
            obj.gev = zeros(1, obj.Ne);
            mznames = fieldnames(obj.mzs);

            for i = 1:obj.Nmzs
                obj.gev(obj.ezi(:, obj.mzs.(mznames{i}).zi)) = ...
                    obj.mzs.(mznames{i}).ev;
            end

            % Jacobiab transpose matrix of each elements
            obj.JT = [acoefs(2, :); acoefs(3, :); acoefs(4, :); ...
                bcoefs(2, :); bcoefs(3, :); bcoefs(4, :); ...
                ccoefs(2, :); ccoefs(3, :); ccoefs(4, :); ];
            [i, j] = getij(3, xNe);
            obj.JT = sparse(i, j, obj.JT(:));
            % getting specefied edata
            edata = getedata(obj.etype);

            switch obj.etype
                case 'TTL4'
                    % evaluation of grad phi i
                    gphi = obj.JT \ repmat(edata.GG, xNe, 1);
                    % inserting in gradx and grady phi i
                    obj.gphix = (gphi(1:3:end, :))';
                    obj.gphiy = (gphi(2:3:end, :))';
                    obj.gphiz = (gphi(3:3:end, :))';
                    % evaluation of Ke
                    obj.Ke.x = zeros(10, xNe);
                    obj.Ke.y = zeros(10, xNe);
                    obj.Ke.z = zeros(10, xNe);
                    temp = 0;

                    for i = 1:4

                        for j = 1:i
                            temp = temp + 1;
                            obj.Ke.x(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :);
                            obj.Ke.y(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :);
                            obj.Ke.z(temp, :) = obj.gphiz(i, :) .* obj.gphiz(j, :);
                        end

                    end

                    obj.Ke.yx = zeros(16, xNe);
                    obj.Ke.zx = zeros(16, xNe);
                    obj.Ke.zy = zeros(16, xNe);

                    for i = 1:4

                        for j = 1:4
                            temp = temp + 1;
                            obj.Ke.yx(temp, :) = obj.gphiy(i, :) .* obj.gphix(j, :);
                            obj.Ke.zx(temp, :) = obj.gphiz(i, :) .* obj.gphix(j, :);
                            obj.Ke.zy(temp, :) = obj.gphiz(i, :) .* obj.gphiy(j, :);
                        end

                    end

                    % multiplying by triangle areas
                    obj.Ke.x = obj.Ke.x * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.y = obj.Ke.y * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.z = obj.Ke.z * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.yx = obj.Ke.yx * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.zx = obj.Ke.zx * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.zy = obj.Ke.zy * sparse(1:xNe, 1:xNe, obj.gev);
                    % elemental force matrix
                    obj.Fe = repmat((obj.gev / 4), 4, 1);
                otherwise
                    error('Element type does not defined.');
            end

            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
        end

        function set.etype(obj, etype)

            if ~ ischar(etype)
                error('Element type must be char.');
            end

            etype = upper(etype);

            if ~ ismember(etype, {'TTL4', 'TTL8'})
                error('Element type does not defined.');
            end

            obj.etype = etype;
        end

        %% Topological Functions
        % setting needed data
        function obj = setdata(obj)
            % first facet of each triangle
            f1 = obj.cl(:, [1, 2, 3]);
            % second facet of each triangle
            f2 = obj.cl(:, [2, 4, 3]);
            % third facet of each triangle
            f3 = obj.cl(:, [3, 4, 1]);
            % forth facet of each triangle
            f4 = obj.cl(:, [1, 4, 2]);
            % sorting for lower index
            [f1, s1] = sort(f1, 2);
            [f2, s2] = sort(f2, 2);
            [f3, s3] = sort(f3, 2);
            [f4, s4] = sort(f4, 2);
            % specefying changed facet index
            s1 = ((s1(:, 1) == 1) & (s1(:, 2) == 3)) | ...
                ((s1(:, 1) == 3) & (s1(:, 2) == 2)) | ...
                ((s1(:, 1) == 2) & (s1(:, 2) == 1));
            s2 = ((s2(:, 1) == 1) & (s2(:, 2) == 3)) | ...
                ((s2(:, 1) == 3) & (s2(:, 2) == 2)) | ...
                ((s2(:, 1) == 2) & (s2(:, 2) == 1));
            s3 = ((s3(:, 1) == 1) & (s3(:, 2) == 3)) | ...
                ((s3(:, 1) == 3) & (s3(:, 2) == 2)) | ...
                ((s3(:, 1) == 2) & (s3(:, 2) == 1));
            s4 = ((s4(:, 1) == 1) & (s4(:, 2) == 3)) | ...
                ((s4(:, 1) == 3) & (s4(:, 2) == 2)) | ...
                ((s4(:, 1) == 2) & (s4(:, 2) == 1));
            % unification of facets
            [obj.facets, ~, ic] = unique([f1; f2; f3; f4], 'rows');
            % getting number of elements
            ne = obj.Ne;
            % getting index of facets corresponding to each elements
            f1 = ic(1:ne);
            f2 = ic(1 + ne:2 * ne);
            f3 = ic(1 + 2 * ne:3 * ne);
            f4 = ic(1 + 3 * ne:4 * ne);
            % specefying boundary facets
            obj.bfacets = sparse([f1, f2, f3, f4], ones(4 * ne, 1), ones(4 * ne, 1));
            obj.bfacets = full(obj.bfacets == 1);
            % specefying trace direction
            f1(s1) = -f1(s1);
            f2(s2) = -f2(s2);
            f3(s3) = -f3(s3);
            f4(s4) = -f4(s4);
            % element matrix
            obj.elements(:, 1:4) = [f1, f2, f3, f4];
            % edge element
            obj.facets = [obj.facets, zeros(size(obj.facets, 1), 2)];
            ttmdbc_evalfe(obj.facets, obj.elements);
            obj.facetNamedSelections.('none') = find(obj.bfacets);
        end

        %% mesh visiualization
        function varargout = showwf(obj)
            obj.ggmesh;

            f = GraphicWindow();
            f.Name = 'Wire Frame Mesh';
            h = guihandles(f);
            delete(h.bg);



            index = obj.facets(:, 4) ~= obj.facets(:, 5);
            patch('Faces', obj.facets(index, 1:3), 'Vertices', ...
                obj.nodes, 'FaceColor', ...
                'c', 'EdgeColor', 'w', ...
                'FaceAlpha', 0.5, 'parent', h.va);
            set(f, 'HandleVisibility', 'off', 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end


        end

        function varargout = showfb(obj)
            obj.ggmesh;



            f = GraphicWindow();
            %             f.Name = ['[Number of Mesh Zones: ', num2str(numel(mzNames)), ']'];
            h = guihandles(f);

            patch('Faces', obj.facets(obj.bfacets, 1:3), 'Vertices', obj.nodes, ...
                'FaceColor', 'b', 'EdgeColor', 'k', ...
                'FaceAlpha', 0.1, 'parent', h.va);
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function f = showm(obj)
            f = GraphicWindow();
            f.Name = 'Global Mesh';
            h = guihandles(f);
            mzNames = fieldnames(obj.mzs);

            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                mzptr.setdata;

                if isequal(mzptr.color, 'none')
                    edgeColor = 'none';
                else
                    edgeColor = [0.1, 0.1, 0.1];
                end

                patch(h.va, 'Faces', mzptr.facets(mzptr.bfacets, :), ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    'c', 'EdgeColor', edgeColor, ...
                    'FaceAlpha', 0.8, 'tag', mzNames{i}, 'ButtonDownFcn', @selectPatchCallback);
            end

            set(f, 'Visible', 'on');
        end

        function showmzs(obj)
            mzNames = fieldnames(obj.mzs);
            f = GraphicWindow();
            f.Name = ['[Number of Mesh Zones: ', num2str(numel(mzNames)), ']'];
            h = guihandles(f);

            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                mzptr.setdata;

                if isequal(mzptr.color, 'none')
                    edgeColor = 'none';
                else
                    edgeColor = [0.1, 0.1, 0.1];
                end

                patch(h.va, 'Faces', mzptr.facets(mzptr.bfacets, :), ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    mzptr.color, 'EdgeColor', edgeColor, ...
                    'FaceAlpha', mzptr.transparency, 'tag', mzNames{i});
            end

            set(f, 'Visible', 'on');
        end

        function showmz(obj, mzName)
            mzNames = fieldnames(obj.mzs);
            a = setFigure(['[Number of Mesh Zones: ', num2str(numel(mzNames)), ']']);

            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                mzptr.setdata;
                patch('Faces', mzptr.facets(mzptr.bfacets, :), ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    'none', 'EdgeColor', 'k', ...
                    'Parent', a, 'FaceAlpha', 0.1, 'edgeAlpha', 0.1);
            end

            mzptr = obj.mzs.(mzName);
            mzptr.setdata;
            patch('Faces', mzptr.facets(mzptr.bfacets, :), ...
                'Vertices', mzptr.nodes, 'FaceColor', ...
                mzptr.color, 'EdgeColor', [0.1, 0.1, 0.1], ...
                'Parent', a, 'FaceAlpha', 1);
            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end

        function showg(obj)
            mzNames = fieldnames(obj.mzs);
            f = GraphicWindow();
            f.Name = ['[Number of Mesh Zones: ', num2str(numel(mzNames)), ']'];
            h = guihandles(f);

            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                mzptr.setdata;

                if isequal(mzptr.color, 'none')
                    edgeColor = 'none';
                else
                    edgeColor = [0.1, 0.1, 0.1];
                end
                
                patch(h.va, 'Faces', mzptr.facets(mzptr.bfacets, :), ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    mzptr.color, 'EdgeColor', mzptr.color, ...
                    'FaceAlpha', mzptr.transparency, 'tag', mzNames{i});
            end

            set(f, 'Visible', 'on');
        end

        function showmzss(obj)
            mznames = fieldnames(obj.mzs);
            ah = setFigure(['Number of Mesh Zones: ', num2str(numel(mznames))]);

            for i = 1:numel(mznames)
                if strcmpi(obj.mzs.(mznames{i}).material, 'air'), continue; end
                obj.mzs.(mznames{i}).setdata;
                patch('Faces', obj.mzs.(mznames{i}).facets(obj.mzs.(mznames{i}).bfacets, :), ...
                    'Vertices', obj.mzs.(mznames{i}).nodes, 'FaceColor', ...
                    obj.mzs.(mznames{i}).color, 'EdgeColor', 'k', ...
                    'Parent', ah, 'FaceAlpha', 1);
            end

            set(gca, 'HandleVisibility', 'off');
        end

        function showElement(obj, eindex)
            index = [obj.e(eindex, 1:4), obj.e(eindex, 6:end)];
            hold all

            for i = 1:length(index)
                plot3(obj.nodes(index(i), 1), obj.nodes(index(i), 2), obj.nodes(index(i), 3), ...
                    'o', 'color', 'y', 'linewidth', 1, 'MarkerFaceColor', 'y');
                text(obj.nodes(index(i), 1), obj.nodes(index(i), 2), obj.nodes(index(i), 3), num2str(i));
            end

            axis off equal

            view([1, 1, 1])
        end

        function varargout = shownfs(obj)

            f = GraphicWindow();
            f.Name = '[Named Facets]';
            h = guihandles(f);

            nfs = fieldnames(obj.facetNamedSelections);

            for i = 1:numel(nfs)
                tmp = rand(1, 3);
                patch('Faces', obj.facets(obj.facetNamedSelections.(nfs{i}), 1:3), 'Vertices', ...
                    obj.nodes, 'FaceColor', ...
                    tmp, 'EdgeColor', 'k', ...
                    'FaceAlpha', 1, 'parent', h.va);
            end

            legend(h.va, nfs);
            set(f, 'HandleVisibility', 'off', 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = shownf(obj, nfName)

            nfName = obj.checkFacetNamedSelectionExistence(nfName);

            f = GraphicWindow();
            f.Name = ['[Named Facet: ', nfName, ']'];
            h = guihandles(f);

            tmp = rand(1, 3);
            patch('Faces', obj.facets(obj.facetNamedSelections.(nfName), 1:3), 'Vertices', ...
                obj.nodes, 'FaceColor', ...
                tmp, 'EdgeColor', 'k', ...
                'FaceAlpha', 1, 'parent', h.va);
            legend(h.va, nfName);

            set(f, 'HandleVisibility', 'off', 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

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

        function shiftMeshZones(obj, mzNames, varargin)
            for mzName = mzNames
                obj.shiftMeshZone(mzName, [shiftX, shiftY]);
            end
        end

        function shmz(varargin)
            shiftMeshZone(varargin{:});
        end

        %% Tools: some operations on mesh zones
        function joinMeshZones(obj, newMeshZoneName, varargin)

            % get mesh zone names string list
            mzNames = emdlab_flib_varargin2StringList(varargin{:});
            Nmzs = numel(mzNames);

            if Nmzs < 1
                error('You have to specify one mesh zone at least.');
            end

            % find total number of mesh zones need to be joined            
            for i = 1:Nmzs
                mzNames(i) = obj.checkMeshZoneExistence(mzNames(i));
            end            

            % store number of nodes and mesh zones of each element
            Nn_tmp = zeros(1, Nmzs);
            Ne_tmp = zeros(1, Nmzs);

            for i = 1:Nmzs
                Nn_tmp(i) = obj.mzs.(mzNames(i)).Nn;
                Ne_tmp(i) = obj.mzs.(mzNames(i)).Ne;
            end

            n_nmz = zeros(sum(Nn_tmp), 3);
            e_nmz = zeros(sum(Ne_tmp), 4);
            n_tmp = 0;
            e_tmp = 0;

            for i = 1:Nmzs
                n_nmz(1 + n_tmp:n_tmp + Nn_tmp(i), :) = obj.mzs.(mzNames(i)).nodes;
                e_nmz(1 + e_tmp:e_tmp + Ne_tmp(i), :) = obj.mzs.(mzNames(i)).cl + n_tmp;
                n_tmp = n_tmp + Nn_tmp(i);
                e_tmp = e_tmp + Ne_tmp(i);
            end

            % jointing mzs
            [n_nmz, ~, ic] = uniquetol(n_nmz, obj.gleps, 'ByRows', true);
            e_nmz = ic(e_nmz);

            % get material and color of first mesh zone
            material = obj.mzs.(mzNames(1)).material;
            color = obj.mzs.(mzNames(1)).color; 
            transparency = obj.mzs.(mzNames(1)).transparency;

            % removing old mesh zones
            for i = 1:Nmzs
                obj.mzs = rmfield(obj.mzs, mzNames(i));
            end

            % check non existance of new mesh zone name
            newMeshZoneName = obj.checkMeshZoneNonExistence(newMeshZoneName);

            % adding new mz
            obj.mzs.(newMeshZoneName) = emdlab_m3d_ttmz(e_nmz, n_nmz);
            obj.mzs.(newMeshZoneName).material = material;
            obj.mzs.(newMeshZoneName).color = color;  
            obj.mzs.(newMeshZoneName).transparency = transparency;

        end

        function jmzs(varargin)
            joinMeshZones(varargin{:});
        end

        function getQuality(obj)
            obj.ggmesh;
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

        function gq(varargin)
            getQuality(varargin{:});
        end

        %% Named Selections
        % facet
        function name = checkFacetNamedSelectionExistence(obj, name)
            name = rmspaces(name);

            if ~ isfield(obj.facetNamedSelections, name)
                error('Specified facet named selection does not exist.');
            end

        end

        function name = checkFacetNamedSelectionNonExistence(obj, name)
            name = rmspaces(name);

            if isfield(obj.facetNamedSelections, name)
                error('Specified facet named selection already exist.');
            end

        end

        function addFacetNamedSelection(obj, name, indices)

            if ~ isvector(indices)
                error('indices must be a column vector.');
            end

            name = obj.checkFacetNamedSelectionNonExistence(name);
            obj.facetNamedSelections.(name) = indices;
            obj.facetNamedSelections.('none') = setdiff(...
                obj.facetNamedSelections.('none'), ...
                obj.facetNamedSelections.(name));
        end

        %% Index Finding
        function y = getfbf(obj)
            obj.ggmesh;
            y = find(obj.bfacets);
        end

        function y = getfbn(obj)
            y = obj.getfbf;
            y = obj.facets(y, 1:3);
            y = unique(y(:));
        end

        function y = getfbfiop(obj, varargin)
            y = ttmdbc_getfbfiop(obj.facets, obj.nodes, obj.facetNamedSelections.none, varargin{:});
            y = y(y ~= 0);
        end

        function y = getfbfiohp(obj, varargin)
            y = ttmdbc_getfbfiohp(obj.facets, obj.nodes, obj.facetNamedSelections.none, varargin{:});
            y = y(y ~= 0);
        end

        function y = getbfioic(obj, p0, u, r)
            y = ttmdbc_getbfioic(obj.facets, obj.nodes, obj.facetNamedSelections.none, p0, u, r);
            y = y(y ~= 0);
        end

        function y = getbfioc(obj, p0, u, r, zb, zu)
            y = ttmdbc_getbfioc(obj.facets, obj.nodes, obj.facetNamedSelections.none, p0, u, r, zb, zu);
            y = y(y ~= 0);
        end

        function y = getbfiohp(obj, p0, u1, u2)
            n = n / norm(n);
            y = ttmdbc_getbfiop(obj.facets, obj.nodes, obj.facetNamedSelections.none, p0, n);
            y = y(y ~= 0);
        end

        function [km, ks] = splitShift(obj, k, vec)
            Nk = length(k);

            if rem(Nk, 2) ~= 0
                error('number of input indices must be even.');
            end

            tp = obj.nodes(k, :);
            km = zeros(Nk / 2, 1);
            ks = zeros(Nk / 2, 1);
            temp = 1;

            while ~ isempty(tp)
                sp = tp(1, :);
                tp = tp(2:end, :);
                km(temp) = k(1);
                k = k(2:end);
                sp = sp + vec;
                index = find(sqrt(sum([tp(:, 1) - sp(1), tp(:, 2) - sp(2), tp(:, 3) - sp(3)].^2, 2)) < obj.geps);

                if index
                    ks(temp) = k(index);
                    tp = tp([1:index - 1, index + 1:end], :);
                    k = k([1:index - 1, index + 1:end]);
                    temp = temp + 1;
                    continue
                end

                sp = sp - 2 * vec;
                index = find(sqrt(sum([tp(:, 1) - sp(1), tp(:, 2) - sp(2), tp(:, 3) - sp(3)].^2, 2)) < obj.geps);

                if index
                    ks(temp) = km(temp);
                    km(temp) = k(index);
                    tp = tp([1:index - 1, index + 1:end], :);
                    k = k([1:index - 1, index + 1:end]);
                    temp = temp + 1;
                    continue
                end

                error('These set of points does not form a set of peridic points.');
            end

        end

        function y = getnIndexOnPlane(obj, p0, n)
            y = find(abs(obj.nodes * n' - p0 * n') < obj.gleps);
        end

        function y = getfb(obj)
            y = obj.facets(obj.bfacets, 1:3);
        end

        function y = geteindex(obj, k)

            y = ismember(obj.facets(:, 1), k) & ismember(obj.facets(:, 2), k) & ismember(obj.facets(:, 3), k);
            y = obj.facets(y, 1:3);

        end

        %% auxiliary functions
        % operation sequence: copy mirror -> copy rotate
        % there are two outputs
        % output is a list containing the name of new mesh zones
        function varargout = aux_cmcr(obj, mzName, xyz0mp, dxyz0mp, xyz0ra, dxyz0ra, Ncopy, rotAngle)

            if nargin < 8
                rotAngle = 2*pi/Ncopy;
            end

            % change the name of current mesh zone
            obj.changeMeshZoneName(mzName, [mzName, '11']);

            % mirror mesh zone with respect to x axis
            obj.cmmz([mzName, '21'], [mzName, '11'], xyz0mp, dxyz0mp);

            for i = 2:Ncopy
                obj.crmz([mzName, '1', num2str(i)], [mzName, '1', num2str(i-1)], rotAngle, xyz0ra, dxyz0ra);
                obj.crmz([mzName, '2', num2str(i)], [mzName, '2', num2str(i-1)], rotAngle, xyz0ra, dxyz0ra);
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

        function aux_cmxyjcmzxjcrzj(obj, mzName, varargin)
            mzName = char(mzName);
            nmzName = [mzName, 'copy_z'];
            obj.copyMirrorMeshZone(nmzName, mzName, [0,0,1]);
            obj.joinMeshZones(mzName, mzName, nmzName);
            newNames = obj.aux_cmcr(mzName, [0,0,0], [0,1,0], [0,0,0], [0,0,1], varargin{:});
            obj.joinMeshZones(mzName, newNames);
        end

        function varargout = aux_cmxyjcmzxcrz(obj, mzName, varargin)
            mzName = char(mzName);
            nmzName = [mzName, 'copy_z'];
            obj.copyMirrorMeshZone(nmzName, mzName, [0,0,1]);
            obj.joinMeshZones(mzName, mzName, nmzName);
            newNames = obj.aux_cmcr(mzName, [0,0,0], [0,1,0], [0,0,0], [0,0,1], varargin{:});

            if nargout == 1, varargout{1} = newNames;
            elseif nargout > 1, error('Too many output arguments.'); end
        end

        %% flag functions
        function makeFalse_isGlobalMeshGenerated(obj)
            obj.isGlobalMeshGenerated = false;
            obj.isJITEvaluated = false;
            obj.isKeMeFe_TTL4_Evaluated = false;
        end

    end

end

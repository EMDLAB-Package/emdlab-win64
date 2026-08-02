% EMDLAB: Electrical Machines Design Laboratory
% face class for 2d geometries
% a face consists several loops, the first loop is outer loop and the rests
% are inner loops

classdef emdlab_g2d_face < handle

    properties

        % a cell to store loops
        loops (1,:) emdlab_g2d_loop;

        % the object tag
        tag (1,:) char = '';

        % constraint mesh points
        meshPoints (:,2) double;

        % color
        color (1,3) double = rand(1,3);

    end

    properties (Dependent = true)

        % number of loops
        Nloops (1,1) double;

    end

    methods

        function obj = emdlab_g2d_face()
        end

        function addLoop(obj, newLoop)
            obj.loops(end+1) = newLoop;
        end

        function y = get.Nloops(obj)
            y = numel(obj.loops);
        end

        function [f, v] = getFacesVertices(obj)

            Nl = numel(obj.loops);
            f = cell(Nl,1);
            v = cell(Nl,1);
            for i = 1:Nl
                [v{i}, f{i}] = obj.loops(i).getMeshNodes;
                if i > 1
                    f{i} = fliplr(f{i});
                end
            end
            Index = 0;
            for i = 2:Nl
                Index = Index + size(v{i-1},1);
                f{i} = f{i} + Index;
            end

            f = cell2mat(f);
            v = cell2mat(v);
            %             [v, ~, ic] = uniquetol(v, 1e-6, 'ByRows', true);
            %             f = ic(f);

        end

        function [f, v] = getFacesVerticesMinimal(obj)

            Nl = numel(obj.loops);
            f = cell(Nl,1);
            v = cell(Nl,1);
            for i = 1:Nl
                [v{i}, f{i}] = obj.loops(i).getMeshNodesMinimal;
                if i > 1
                    f{i} = fliplr(f{i});
                end
            end
            Index = 0;
            for i = 2:Nl
                Index = Index + size(v{i-1},1);
                f{i} = f{i} + Index;
            end

            f = cell2mat(f);
            v = cell2mat(v);

        end

        function m = getMesh(obj, meshGenerator)

            % get faces and vertices
            [f, v] = obj.getFacesVertices;
            v = [v;obj.meshPoints];

            % call mesh generator
            switch meshGenerator
                case 'mm'
                    m = emdlab_m2d_mm(f, v);
                case 'mg0'
                    [f1, v1] = obj.getFacesVerticesMinimal;
                    m = emdlab_m2d_mg0(f, v, f1, v1);
                 case 'mg3'
                    [f1, v1] = obj.getFacesVerticesMinimal;
                    m = emdlab_m2d_mg3(f, v, f1, v1);
            end

        end

        function m = getMeshNew(obj)

            % get faces and vertices
            [f, v] = obj.getFacesVertices;
            m = emdlab_m2d_mg0(f, v, obj.getPolyshape);

        end

        function p = getPolyshape(obj)

            p = polyshape(obj.loops(1).getMeshNodes);
            for i = 2:numel(obj.loops)
                p = p.addboundary(obj.loops(i).getMeshNodes);
            end

        end
    
        function addMeshPoints(obj,p)

            obj.meshPoints = [obj.meshPoints;p];

        end

    end

end

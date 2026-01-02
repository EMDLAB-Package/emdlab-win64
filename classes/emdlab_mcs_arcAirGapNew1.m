% EMDLAB: Electrical Machines Design Laboratory
% moving contacts
% arc air gap: two concentric arcs: pure cylindrical motion

classdef emdlab_mcs_arcAirGapNew1 < handle & emdlab_g2d_constants & matlab.mixin.SetGet
    
    properties (SetAccess = private)
        
        % mesh
        m (1,1) emdlab_m2d_tmz = emdlab_m2d_tmz([], []);

        % inner points
        ips (:,2) double;

        % outer points
        ops (:,2) double;

        % inner point angles
        ipas (:,1) double;

        % outer point angles
        opas (:,1) double;

        % arc center
        center (1,2) double;

        % minimum mesh length
        h0 (1,1) double;

        % number of layers for mesh
        Nlayer (1,1) double {mustBeInteger, mustBeNonnegative} = 2;

        % inner circle radius
        rin (1,1) double;

        % outer circle radius
        rout (1,1) double;

        % number of inner points
        Nips (1,1) double;

        % number of outer points
        Nops (1,1) double;

        % arc angle
        angle (1,1) double;

        % moving boundary
        movingBoundary = 'inner';

        imesh;
        omesh;
        
    end
    
    methods
        
        function obj = emdlab_mcs_arcAirGapNew1(center, ips, ops, Nlayer, movingBoundary)
            
            if nargin < 5, movingBoundary = 'inner'; end
            obj.center = center;
            obj.ips = ips;
            obj.ops = ops;
            
            obj.Nips = size(ips, 1);
            obj.Nops = size(ops, 1);
            
            obj.rin = obj.get_rin;
            obj.rout = obj.get_rout;
            
            obj.Nlayer = Nlayer;
            obj.movingBoundary = movingBoundary;            

            % chekers: inner points must be insode outer points
            obj.sortPoints;

            % evaluation of h0
            l1 = sqrt(sum((obj.ips(1:end - 1, :) - obj.ips(2:end, :)).^2, 2));
            l2 = sqrt(sum((obj.ops(1:end - 1, :) - obj.ops(2:end, :)).^2, 2));
            obj.h0 = mean([l1; l2]);
            obj.updateMesh;
            
        end
        
        function set.Nlayer(obj, newValue)
            
            if newValue < 1
                throw(MException('', 'Number of layers must be higher than 1.'));
            end
            
            obj.Nlayer = newValue;

        end
        
        function y = get_rin(obj)
            
            % calculation of inner circle radius
            tmp = [obj.ips(:, 1) - obj.center(1), obj.ips(:, 2) - obj.center(2)];
            tmp = sqrt(sum(tmp.^2, 2));
            y = mean(tmp);
            
            if sum(abs(tmp - y)) > obj.Nips * obj.gleps
                error('Inner points do not form an arc.');
            end
            
        end
        
        function y = get_rout(obj)
            
            % calculation of inner circle radius
            tmp = [obj.ops(:, 1) - obj.center(1), obj.ops(:, 2) - obj.center(2)];
            tmp = sqrt(sum(tmp.^2, 2));
            y = mean(tmp);
            
            if sum(abs(tmp - y)) > obj.Nops * obj.gleps
                error('Outer points do not form an arc.');
            end
            
        end
        
        function sortPoints(obj)
            
            % inner points
            obj.ipas = atan_02pi([obj.ips(:, 1) - obj.center(1), obj.ips(:, 2) - obj.center(2)]);
            [~, index] = emdlab_g2d_sortPointsCCW(obj.ips, obj.center);
            obj.ips = obj.ips(index, :);
            obj.ipas = obj.ipas(index, :);
            % inner arc angle: as a default it assigned to arc angle
            if obj.ipas(1) >= 0 && obj.ipas(1) < pi
                obj.angle = obj.ipas(end) - obj.ipas(1);
            else
                obj.angle = 2 * pi - obj.ipas(1) + obj.ipas(end);
            end
            
            % outer points
            obj.opas = atan_02pi([obj.ops(:, 1) - obj.center(1), obj.ops(:, 2) - obj.center(2)]);
            [~, index] = emdlab_g2d_sortPointsCCW(obj.ops, obj.center);
            obj.ops = obj.ops(index, :);
            obj.opas = obj.opas(index, :);
            % outer arc angle
            if obj.opas(1) >= 0 && obj.opas(1) < pi
                tmp = obj.opas(end) - obj.opas(1);
            else
                tmp = 2 * pi - obj.opas(1) + obj.opas(end);
            end
            
            % arc angles must be the same with a geometry tolerance
            if abs(obj.angle - tmp) > obj.gaeps
                error('Arc angles must be the same.');
            end
            
        end
        
        function updateMesh(obj)
            
            % finding priority
            outerProduct = det([obj.ops(1,:);obj.ips(1,:)])/obj.rin/obj.rout;
            innerProduct = sum(obj.ops(1,:).*obj.ips(1,:))/obj.rin/obj.rout;
            
            % relative angle between starting point of inner arc respect to outer arc
            if innerProduct > 0
                alpha = abs(asin(outerProduct));
            else
                alpha = pi - abs(asin(outerProduct)); 
            end
            
            % mid radius
            rmid = (obj.rin + obj.rout) / 2;
            
            if alpha < obj.gaeps

                % number of points on inner arc layers
                Ntmp = ceil(rmid * obj.angle / obj.h0);

                rtmp = linspace(obj.rin,obj.rout,obj.Nlayer+1);
                rtmp = rtmp(2:end-1);
                
                if obj.ipas(1)<pi
                    a1 = obj.ipas(1);
                else
                    a1 = -2*pi + obj.ipas(1);
                end
                
                t = linspace(a1, obj.ipas(end), Ntmp)';

                p = obj.ips;                
                for i = 1:obj.Nlayer-1
                    p = [p;rtmp(i) * [cos(t), sin(t)]];
                end
                p = [p;obj.ops];
                
                % indicies of critical points
                n = [obj.Nips,Ntmp*ones(1,obj.Nlayer-1),obj.Nops];
                n1 = cumsum([1,n(1:end-1)]);
                n2 = cumsum(n);
                ices = zeros([],2);

                for i = 2:obj.Nlayer
                    newIndex = [n1(i):n2(i)-1;n1(i)+1:n2(i)]';
                    ices = [ices;newIndex];
                end
                
                newIndex = [n1(end):n2(end)-1;n1(end)+1:n2(end)]';
                f = newIndex;
                for i = 1:obj.Nlayer
                    f = [f;[n2(end-i+1),n2(end-i)]];
                end
                newIndex = [n1(1):n2(1)-1;n1(1)+1:n2(1)]';
                f = [f;fliplr(newIndex)];
                for i = 1:obj.Nlayer
                    f = [f;[n1(i),n1(i+1)]];
                end

                obj.m = emdlab_m2d_mm(f, p, ices);

%                 tmp = obj.m.getCenterOfElements;
%                 tmp = tmp - obj.center;
%                 tmp = sqrt(sum(tmp.^2,2));
%                 tmp = tmp<(obj.rin+obj.rout)/2;
% 
%                 obj.imesh = triangulation(obj.m.cl(tmp,:),obj.m.nodes);
%                 obj.omesh = triangulation(obj.m.cl(~tmp,:),obj.m.nodes);
% 
%                 obj.imesh = emdlab_
                return;
                
            else
                
                % number of points on non-overlap arcs
                NtmpSide = ceil(rmid * alpha / obj.h0);
                if NtmpSide < 2, NtmpSide = 2; end
                
                % number of points on overlap arc
                NtmpIn = ceil(rmid * (obj.angle - alpha) / obj.h0);
                if NtmpIn < 2, NtmpIn = 2; end
                
                if outerProduct > 0
                    
                    if obj.opas(1)<pi
                        oa1 = obj.opas(1);
                    else
                        oa1 = -2*pi + obj.opas(1);
                    end
                    
                    if obj.ipas(1)<pi
                        ia1 = obj.ipas(1);
                    else
                        ia1 = -2*pi + obj.ipas(1);
                    end
                
                    % calculation of coordinates of points on mid arc
                    t = linspace(oa1, ia1, NtmpSide);
                    t = [t(1:end - 1), linspace(ia1, obj.opas(end), NtmpIn)];
                    t = [t(1:end - 1), linspace(obj.opas(end), obj.ipas(end), NtmpSide)];
                    p = rmid * [cos(t); sin(t)]';
                    
                    % indicies of critical points
                    index1 = 1;
                    index2 = obj.Nops;
                    index3 = obj.Nips + obj.Nops + NtmpSide + NtmpIn - 1;
                    index4 = obj.Nips + obj.Nops + 2*NtmpSide + NtmpIn - 2;
                    index5 = obj.Nips + obj.Nops;
                    index6 = obj.Nops + 1;
                    index7 = obj.Nips + obj.Nops + NtmpSide;
                    index8 = obj.Nips + obj.Nops + 1;
                    
                    % geometry faces
                    f = [index1:index2-1; index1+1:index2]';
                    f = [f; [index2,index3]];
                    f = [f; [index3:index4-1; index3+1:index4]'];
                    f = [f; [index4,index5]];
                    f = [f; [index5:-1:index6+1; index5-1:-1:index6]'];
                    f = [f; [index6,index7]];
                    f = [f; [index7:-1:index8+1; index7-1:-1:index8]'];
                    f = [f; [index8,index1]];
                    
                    % geometry constraint edges
                    ices = [index7:index3 - 1; index7 + 1:index3]';
                    
                else
                    
                    if obj.opas(1)<pi
                        oa1 = obj.opas(1);
                    else
                        oa1 = -2*pi + obj.opas(1);
                    end
                    
                    if obj.ipas(1)<pi
                        ia1 = obj.ipas(1);
                    else
                        ia1 = -2*pi + obj.ipas(1);
                    end
                
                    % calculation of coordinates of points on mid arc
                    t = linspace(ia1, oa1, NtmpSide);
                    t = [t(1:end - 1), linspace(oa1, obj.ipas(end), NtmpIn)];
                    t = [t(1:end - 1), linspace(obj.ipas(end), obj.opas(end), NtmpSide)];
                    p = rmid * [cos(t); sin(t)]';
                    
                    % indicies of critical points
                    index1 = 1;
                    index2 = obj.Nops;
                    index3 = obj.Nips + obj.Nops + 2*NtmpSide + NtmpIn - 2;
                    index4 = obj.Nips + obj.Nops + NtmpSide + NtmpIn - 1;
                    index5 = obj.Nips + obj.Nops;
                    index6 = obj.Nops + 1;
                    index7 = obj.Nips + obj.Nops + 1;
                    index8 = obj.Nips + obj.Nops + NtmpSide;
                    
                    % geometry faces
                    f = [index1:index2-1; index1+1:index2]';
                    f = [f; [index2,index3]];
                    f = [f; [index3:-1:index4+1; index3-1:-1:index4]'];
                    f = [f; [index4,index5]];
                    f = [f; [index5:-1:index6+1; index5-1:-1:index6]'];
                    f = [f; [index6,index7]];
                    f = [f; [index7:index8-1; index7+1:index8]'];
                    f = [f; [index8,index1]];                    
                    
                    % geometry constraint edges
                    ices = [index8:index4 - 1; index8 + 1:index4]';
                    
                end
                
            end
            
            obj.m = emdlab_m2d_mm(f, [obj.ops; obj.ips; p], ices);
            
        end
        
        function obj = rotateInner(obj, alpha)
            
            alpha = rem(alpha, obj.angle);
            if abs(alpha) < obj.gaeps, return; end
            obj.ipas = obj.ipas + alpha;
            obj.ips = ext_protate2(obj.ips, alpha, obj.center);
            obj.updateMesh;
            
        end
        
        function obj = rotateOuter(obj, alpha)
            
            alpha = rem(alpha, obj.angle);
            if abs(alpha) < obj.gaeps, return; end
            obj.opas = obj.opas + alpha;
            obj.ops = ext_protate2(obj.ops, alpha, obj.center);
            obj.updateMesh;
            
        end
        
        function showArcs(obj)
            
            figure; 
            hold on;
            plot(obj.ips(:,1), obj.ips(:,2));
            plot(obj.ops(:,1), obj.ops(:,2));
            
        end
        
        function rotate(obj, varargin)
            if strcmpi(obj.movingBoundary, 'inner')
                obj.rotateInner(varargin{:});
            else
                obj.rotateOuter(varargin{:});
            end
        end

%         function rotateSliding(obj, varargin)
%             if strcmpi(obj.movingBoundary, 'inner')
%                 obj.imesh.
%             else
%                 obj.rotateOuter(varargin{:});
%             end
%         end

    end
    
end

% moving contacts -> line air-gap
% it used for full air-gap modeling of linear electrical machines
% this object is for meshing of a linear band

classdef emdlab_mcs_lineAirGap < handle & emdlab_g2d_constants & matlab.mixin.SetGet
    
    properties (SetAccess = private)
        
         % mesh: triangular mesh zone
        m (1,1);

        % stationary points
        spts (:,2) double;

        % moving points
        mpts (:,2) double;

        % inner point distances
        spds (:,1) double;

        % outer point distances
        mpds (:,1) double;

        % minimum mesh length
        h0 (1,1) double;

        % stationary pivot vector
        svec (1,2) double;

        % moving pivot vector
        mvec (1,2) double;

        % unit pivot vector
        uvec (1,2) double;

        % unit normal vector
        nvec (1,2) double;

        % airgap length
        g (1,1) double;

    end
    properties (Dependent = true)

        % number of inner points
        Nsps (1,1) double;

        % number of outer points
        Nmps (1,1) double;

    end

    methods
        function obj = emdlab_mcs_lineAirGap(sps, mps, varargin)
%             set(obj, varargin{:});
            % chekers: sps and mps must be on line
            obj.spts = sortPointsOnLineX(sps);
            obj.mpts = sortPointsOnLineX(mps);
            obj.svec = obj.spts(end,:) - obj.spts(1,:);
            obj.mvec = obj.mpts(end,:) - obj.mpts(1,:);
            if abs(norm(obj.svec) - norm(obj.mvec)) > 1e-6
                error('Line segments must be the same size');
            end
            obj.uvec = obj.svec/norm(obj.svec);
            obj.nvec = ext_protate2(obj.uvec, pi/2);
            % evaluation of airgap length
            obj.g = sum(obj.nvec .* (obj.mpts(1,:) - obj.spts(1,:)));
            % evaluation of h0
            l1 = sqrt(sum((obj.spts(1:end-1,:) - obj.spts(2:end,:)).^2, 2));
            l2 = sqrt(sum((obj.mpts(1:end-1,:) - obj.mpts(2:end,:)).^2, 2));
            obj.h0 = mean([l1; l2]);
            obj.updateMesh;
        end
        function y = get.Nmps(obj)
            y = size(obj.mpts, 1);
        end
        function y = get.Nsps(obj)
            y = size(obj.spts, 1);
        end
        function updateMesh(obj)
            ces1 = [1:obj.Nsps-1;2:obj.Nsps]';
            ces2 = [1:obj.Nmps-1;2:obj.Nmps]' + obj.Nsps;
            dis = sum(obj.uvec .* (obj.mpts(1,:) - obj.spts(1,:)));
            if abs(dis) < obj.h0
                p1 = obj.spts(1,:) + obj.g/2 * obj.nvec;
                p2 = obj.spts(end,:) + obj.g/2 * obj.nvec;
                Np = ceil(norm(obj.svec)/obj.h0);
                if Np<2, Np = 2; end
                xp = linspace(p1(1), p2(1), Np);
                yp = linspace(p1(2), p2(2), Np);
                ip1 = obj.Nsps + obj.Nmps + 1;
                ip2 = obj.Nsps + obj.Nmps + Np;
                ces = [ces1; fliplr(ces2);...
                    obj.Nsps,ip2; ip2, ip1-1; ip1, 1; obj.Nsps+1, ip1];
                ces3 = [ip1:ip2-1; ip1+1:ip2]';
            elseif dis>0
                p1 = obj.spts(1,:) + obj.g/2 * obj.nvec;
                p2 = obj.mpts(1,:) - obj.g/2 * obj.nvec;
                p3 = obj.spts(end,:) + obj.g/2 * obj.nvec;
                p4 = obj.mpts(end,:) - obj.g/2 * obj.nvec;
                Np_side = ceil(dis/obj.h0);
                if Np_side<2, Np_side = 2; end
                Np_mid = ceil((norm(obj.svec) - dis)/obj.h0);
                if Np_mid<2, Np_mid = 2; end
                xp = linspace(p1(1), p2(1), Np_side);
                yp = linspace(p1(2), p2(2), Np_side);
                xp = [xp(1:end-1), linspace(p2(1), p3(1), Np_mid)];
                yp = [yp(1:end-1), linspace(p2(2), p3(2), Np_mid)];
                xp = [xp(1:end-1), linspace(p3(1), p4(1), Np_side)];
                yp = [yp(1:end-1), linspace(p3(2), p4(2), Np_side)];
                ip1 = obj.Nsps + obj.Nmps + 1;
                ip2 = ip1 - 1 + Np_side;
                ip3 = ip2 - 1 + Np_mid;
                ip4 = ip3 - 1 + Np_side;
                tmp1 = [ip1:ip2-1; ip1+1:ip2]';
                tmp2 = [ip3:ip4-1; ip3+1:ip4]';
                ces = [ces1; fliplr(ces2);...
                    obj.Nsps, ip3; ip4, ip1-1; ip1, 1; obj.Nsps+1, ip2;...
                    fliplr(tmp1); tmp2];
                ces3 = [ip2:ip3-1; ip2+1:ip3]';
            else
                p1 = obj.mpts(1,:) - obj.g/2 * obj.nvec;
                p2 = obj.spts(1,:) + obj.g/2 * obj.nvec;
                p3 = obj.mpts(end,:) - obj.g/2 * obj.nvec;
                p4 = obj.spts(end,:) + obj.g/2 * obj.nvec;
                dis = abs(dis);
                Np_side = ceil(dis/obj.h0);
                if Np_side<2, Np_side = 2; end
                Np_mid = ceil((norm(obj.svec) - dis)/obj.h0);
                if Np_mid<2, Np_mid = 2; end
                xp = linspace(p1(1), p2(1), Np_side);
                yp = linspace(p1(2), p2(2), Np_side);
                xp = [xp(1:end-1), linspace(p2(1), p3(1), Np_mid)];
                yp = [yp(1:end-1), linspace(p2(2), p3(2), Np_mid)];
                xp = [xp(1:end-1), linspace(p3(1), p4(1), Np_side)];
                yp = [yp(1:end-1), linspace(p3(2), p4(2), Np_side)];
                ip1 = obj.Nsps + obj.Nmps + 1;
                ip2 = ip1 - 1 + Np_side;
                ip3 = ip2 - 1 + Np_mid;
                ip4 = ip3 - 1 + Np_side;
                tmp1 = [ip1:ip2-1; ip1+1:ip2]';
                tmp2 = [ip3:ip4-1; ip3+1:ip4]';
                ces = [ces1; fliplr(ces2);...
                    obj.Nsps, ip4; ip3, ip1-1; ip2, 1; obj.Nsps+1, ip1;...
                    tmp1; fliplr(tmp2)];
                ces3 = [ip2:ip3-1; ip2+1:ip3]';
            end
            if obj.g<0
                ces = fliplr(ces);
            end
            obj.m = emdlab_m2d_mm(ces, [obj.spts;obj.mpts;xp(:),yp(:)], ces3);
        end
        function obj = shiftStationaryPoints(obj, value)
            obj.spts = ext_pshift2(obj.spts, value * obj.uvec);
            obj.updateMesh;
        end
        function obj = shiftMovingPoints(obj, value)
            obj.mpts = ext_pshift2(obj.mpts, value * obj.uvec);
            obj.updateMesh;
        end
    end
end

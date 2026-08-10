% moving contacts -> line air-gap
% it used for full air-gap modeling of linear electrical machines
% this object is for meshing of a linear band

classdef emdlab_mcs_lineAirGap1 < handle & emdlab_g2d_constants & matlab.mixin.SetGet
    
    properties (SetAccess = private)
        
        % triangular mesh zone pointer
        m (1,1);

        % stationary points
        spts (:,2) double;

        % moving points
        mpts (:,2) double;

        % mean mesh length
        h0 (1,1) double;

        % unit vector along moving direction
        u (1,2) double;

        % normal vector to u
        n (1,2) double;

        % airgap length
        gap (1,1) double;

        % number of airgap layers
        Nagl (1,1) double;

        % contact length
        L (1,1) double;

    end

    properties (Dependent = true)

        % number of inner points
        Nspts (1,1) double;

        % number of outer points
        Nmpts (1,1) double;

    end

    methods

        function obj = emdlab_mcs_lineAirGap1(spts, mpts, ux, uy, Nagl)

            if nargin < 5, Nagl = 2; end
            tol = 1e-6;
            uMag = norm([ux,uy]);
            ux = ux/uMag;
            uy = uy/uMag;

            % chekers: spts and mpts must be on line
            [~, obj.spts] = emdlab_g2d_arePointsOnLine(spts, spts(1,1), spts(1,2), ux, uy, tol);
            [~, obj.mpts] = emdlab_g2d_arePointsOnLine(mpts, mpts(1,1), mpts(1,2), ux, uy, tol);

            % tangent and normal vectors
            obj.u = [ux,uy];
            obj.n = [-uy,ux];

            % calculation of the airgap length
            obj.gap = obj.n * (obj.mpts(1,:) - obj.spts(1,:))';

            % evaluation of h0
            L1 = vecnorm(diff(obj.spts), 2, 2);
            L2 = vecnorm(diff(obj.mpts), 2, 2);
            obj.h0 = mean([L1; L2]);
            L1 = sum(L1);
            L2 = sum(L2);
            if abs(L1 - L2) > 1e-5
                error('Contact boundaries dont have the same length.');
            end
            obj.L = (L1 + L2)/2;
            obj.Nagl = Nagl;
            obj.mpts = flipud(obj.mpts);
            obj.updateMesh;

        end
        
        function y = get.Nmpts(obj)
            y = size(obj.mpts, 1);
        end

        function y = get.Nspts(obj)
            y = size(obj.spts, 1);
        end

        function updateMesh(obj)

            % constraint edges #1
            ces1 = [1:obj.Nspts-1;2:obj.Nspts]';

            % constraint edges #2
            ces3 = [1:obj.Nmpts-1;2:obj.Nmpts]' + obj.Nspts;

            % distance between first stationary and moving point along u-vec
            dis = sum(obj.u .* (obj.mpts(end,:) - obj.spts(1,:)));

            if abs(dis) < obj.h0

                % number of inner layer points
                Nils = obj.Nagl - 1;

                % number of inner points
                Nipts = max(ceil(obj.L/obj.h0), 2);

                % allocate memory for inner points
                ipts = zeros(Nils * Nipts, 2);

                % side vector
                svec = obj.mpts(end,:) - obj.spts(1,:);
                svec = svec/norm(svec);

                for i = 1:Nils

                    % index of new inner points: new row
                    idx = (1:Nipts) + (i-1) * Nipts;

                    % start and end points of new row
                    p1 = obj.spts(1,:) + (i*obj.gap/obj.Nagl) * svec;
                    p2 = obj.spts(end,:) + (i*obj.gap/obj.Nagl) * svec;

                    % x and y cordinates of new row
                    xr = linspace(p1(1), p2(1), Nipts);
                    yr = linspace(p1(2), p2(2), Nipts);

                    % store new points
                    ipts(idx,:) = [xr;yr]';

                end

                % calculation of ces2
                idx_ip = obj.Nspts + obj.Nmpts + (Nipts:Nipts:(Nipts*Nils));
                idx_ip = [obj.Nspts, idx_ip, obj.Nspts+1];
                idx_ip = [idx_ip(1:end-1);idx_ip(2:end)];
                ces2 = idx_ip';

                % calculation of ces4
                idx_ip = obj.Nspts + obj.Nmpts + (1:Nipts:((Nipts-1)*Nils));
                idx_ip = [1, idx_ip, obj.Nspts+obj.Nmpts];
                idx_ip = fliplr(idx_ip);
                idx_ip = [idx_ip(1:end-1);idx_ip(2:end)];
                ces4 = idx_ip';
                
                % facets and vertices list for mesh generation
                F = [ces1; ces2; ces3; ces4];
                V = [obj.spts; obj.mpts; ipts];

            elseif dis > 0

                % check for overlap
                if dis >= obj.L
                    error('There is no overlap between contact & target.');
                end

                % set number of airgap layers to an even number
                if mod(obj.Nagl,2), obj.Nagl = obj.Nagl + 1; end

                p1 = obj.spts(1,:) + obj.gap/2 * obj.n;
                p2 = obj.mpts(1,:) - obj.gap/2 * obj.n;
                p3 = obj.spts(end,:) + obj.gap/2 * obj.n;
                p4 = obj.mpts(end,:) - obj.gap/2 * obj.n;
                Np_side = ceil(dis/obj.h0);
                if Np_side<2, Np_side = 2; end
                Np_mid = ceil((norm(obj.svec) - dis)/obj.h0);
                if Np_mid<2, Np_mid = 2; end
                xr = linspace(p1(1), p2(1), Np_side);
                yr = linspace(p1(2), p2(2), Np_side);
                xr = [xr(1:end-1), linspace(p2(1), p3(1), Np_mid)];
                yr = [yr(1:end-1), linspace(p2(2), p3(2), Np_mid)];
                xr = [xr(1:end-1), linspace(p3(1), p4(1), Np_side)];
                yr = [yr(1:end-1), linspace(p3(2), p4(2), Np_side)];
                ip1 = obj.Nspts + obj.Nmpts + 1;
                ip2 = ip1 - 1 + Np_side;
                ip3 = ip2 - 1 + Np_mid;
                ip4 = ip3 - 1 + Np_side;
                tmp1 = [ip1:ip2-1; ip1+1:ip2]';
                tmp2 = [ip3:ip4-1; ip3+1:ip4]';
                ces = [ces1; fliplr(ces2);...
                    obj.Nspts, ip3; ip4, ip1-1; ip1, 1; obj.Nspts+1, ip2;...
                    fliplr(tmp1); tmp2];
                ces3 = [ip2:ip3-1; ip2+1:ip3]';

            else

                % check for overlap
                if dis <= -obj.L
                    error('There is no overlap between contact & target.');
                end

                % set number of airgap layers to an even number
                if mod(obj.Nagl,2), obj.Nagl = obj.Nagl + 1; end

                p1 = obj.mpts(1,:) - obj.gap/2 * obj.n;
                p2 = obj.spts(1,:) + obj.gap/2 * obj.n;
                p3 = obj.mpts(end,:) - obj.gap/2 * obj.n;
                p4 = obj.spts(end,:) + obj.gap/2 * obj.n;
                dis = abs(dis);
                Np_side = ceil(dis/obj.h0);
                if Np_side<2, Np_side = 2; end
                Np_mid = ceil((norm(obj.svec) - dis)/obj.h0);
                if Np_mid<2, Np_mid = 2; end
                xr = linspace(p1(1), p2(1), Np_side);
                yr = linspace(p1(2), p2(2), Np_side);
                xr = [xr(1:end-1), linspace(p2(1), p3(1), Np_mid)];
                yr = [yr(1:end-1), linspace(p2(2), p3(2), Np_mid)];
                xr = [xr(1:end-1), linspace(p3(1), p4(1), Np_side)];
                yr = [yr(1:end-1), linspace(p3(2), p4(2), Np_side)];
                ip1 = obj.Nspts + obj.Nmpts + 1;
                ip2 = ip1 - 1 + Np_side;
                ip3 = ip2 - 1 + Np_mid;
                ip4 = ip3 - 1 + Np_side;
                tmp1 = [ip1:ip2-1; ip1+1:ip2]';
                tmp2 = [ip3:ip4-1; ip3+1:ip4]';
                ces = [ces1; fliplr(ces2);...
                    obj.Nspts, ip4; ip3, ip1-1; ip2, 1; obj.Nspts+1, ip1;...
                    tmp1; fliplr(tmp2)];
                ces3 = [ip2:ip3-1; ip2+1:ip3]';

            end

            if obj.gap<0
                F = fliplr(F);
            end

            obj.m = emdlab_m2d_mm(F,V);

        end

        function obj = shiftStationaryPoints(obj, value)
            obj.spts = ext_pshift2(obj.spts, value * obj.u);
            obj.updateMesh;
        end

        function obj = shiftMovingPoints(obj, value)
            obj.mpts = ext_pshift2(obj.mpts, value * obj.u);
            obj.updateMesh;
        end

    end

end

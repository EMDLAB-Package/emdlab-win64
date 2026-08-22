% EMDLAB: Electrical Machines Design Laboratory
% two-dimensional arc class

classdef emdlab_g2d_arc < handle & emdlab_g2d_constants

    properties

        % center, p1, p2 points
        p0 (1,1) emdlab_g2d_point;
        p1 (1,1) emdlab_g2d_point;
        p2 (1,1) emdlab_g2d_point;

        % direction of the arc
        direction (1,1) logical = true;

        Nnodes (1,1) double = 5;
        isSetNnodes (1,1) logical = true;

        maxLength (1,1) double = 1;
        isSetMaxLength (1,1) logical = false;

        maxDegree (1,1) double = deg2rad(20);
        isSetMaxDegree (1,1) logical = false;

        healTol (1,1) double;
        L1 (1,1) double;
        L2 (1,1) double;

    end

    methods

        % constructor
        function obj = emdlab_g2d_arc(p0, p1, p2, direction)

            obj.p0 = p0;
            obj.p1 = p1;
            obj.p2 = p2;

            if nargin<4, direction = true; end
            obj.direction = direction;

            obj.validateArc;

        end

        function validateArc(obj)

            % check for consistency
            p1c = obj.p1 - obj.p0;
            p2c = obj.p2 - obj.p0;
            r1 = p1c.getRadialLength;
            r2 = p2c.getRadialLength;

            if abs(r2-r1) > max(obj.gleps, obj.healTol)

                throw(MException('', 'These points do not form an arc.'));

            end

            if obj.healTol > obj.gleps

                obj.p2.setVector(obj.p0.getVector + r1 * p2c.getUnitVector);

            end

        end

        % setters
        function setNnodes(obj, newNnodes)

            obj.Nnodes = max(newNnodes,3);
            obj.isSetMaxLength = false;
            obj.isSetMaxDegree = false;
            obj.isSetNnodes = true;

        end

        function setMaxLength(obj, newMaxLength)

            R = obj.getRadius;
            if obj.maxDegree > obj.getAngle
                obj.maxDegree = obj.getAngle / 3;
            end
            tmp_length = 2*R*sin(obj.maxDegree*0.5);
            if newMaxLength > tmp_length
                newMaxLength = tmp_length;
            end
            obj.maxLength = newMaxLength;
            obj.p1.meshSize = newMaxLength;
            obj.p2.meshSize = newMaxLength;
            obj.isSetMaxLength = true;
            obj.isSetMaxDegree = false;
            obj.isSetNnodes = false;

        end

        % getters
        function y = getRadius(obj)

            y = norm(obj.p1-obj.p0);

        end

        function y = getAngle(obj)

            p1c = obj.p1 - obj.p0;
            p2c = obj.p2 - obj.p0;
            p2p1 = obj.p2 - obj.p1;
            y = 2*real(asin(norm(p2p1)/2/obj.getRadius));
            if obj.direction
                if p1c.outerProduct(p2c) < 0
                    y = 2*pi - y ;
                end
            else
                if p1c.outerProduct(p2c) > 0
                    y = 2*pi - y ;
                end
            end

        end

        function y = getSignedAngle(obj)
            y = obj.getAngle;
            if ~obj.direction
                y = -y;
            end
        end

        function y = getAngleDegree(obj)

            y = obj.getAngle*180/pi;

        end

        function y = getLength(obj)
            y = obj.getRadius * obj.getAngle;
        end

        %         function y = getCenter(obj)
        %
        %             p = obj.getMeshNodes();
        %             plot(p(:,1),p(:,2));
        %             n = size(p,1);
        %             switch rem(n,2)
        %                 case 0
        %                     p1_tmp = p(n/2,:);
        %                     p2_tmp = p(n/2+1,:);
        %                     y = (p1_tmp +  p2_tmp)/2;
        %                 case 1
        %                     y = p(ceil(n/2),:);
        %             end
        %
        %         end

        function y = getCenter(obj)
            if obj.direction
                y = emdlab_g2d_rotatePoints(obj.p1.getVector,obj.getAngle/2,obj.p0.x,obj.p0.y);
            else
                y = emdlab_g2d_rotatePoints(obj.p1.getVector,-obj.getAngle/2,obj.p0.x,obj.p0.y);
            end
        end

        function nodes = getMeshNodesMinimal(obj)

            nodes = obj.getMeshNodes;

        end

        function [nodes, cl] = getMeshNodes(obj)

            if obj.isSetNnodes

                nodes = zeros(obj.Nnodes, 2);
                stepAngle = obj.getAngle/(obj.Nnodes-1);

                if ~obj.direction
                    stepAngle = -stepAngle;
                end

                c = obj.p0.getVector;
                nodes(1,:) = obj.p1.getVector;
                for i = 2:obj.Nnodes-1
                    nodes(i,:) = ext_protate2(nodes(1,:), (i-1)*stepAngle, c);
                end
                nodes(end,:) = obj.p2.getVector;

            elseif obj.isSetMaxLength

                obj.Nnodes = max(ceil(obj.getRadius*obj.getAngle/obj.maxLength), 2);
                nodes = zeros(obj.Nnodes, 2);
                stepAngle = obj.getAngle/(obj.Nnodes-1);

                if ~obj.direction
                    stepAngle = -stepAngle;
                end

                c = obj.p0.getVector;
                nodes(1,:) = obj.p1.getVector;
                for i = 2:obj.Nnodes-1
                    nodes(i,:) = ext_protate2(nodes(1,:), (i-1)*stepAngle, c);
                end
                nodes(end,:) = obj.p2.getVector;

            end

            if nargout == 2
                % connectivity list
                Nn = size(nodes,1);
                cl = [1:(Nn-1);2:Nn]';
            end

        end

        function y = getInnerPoints(obj)

            y = obj.getMeshPoints;
            y = y(2:end-1,:);

        end

        function [e,v] = getev(obj)
            v = [obj.p1;obj.getip;obj.p2];
            Nv = size(v,1);
            e = [1:Nv-1;2:Nv]';
        end

        function y = getd(obj,p)
            y = min(dsegment(p,[obj.p1;obj.getip;obj.p2]),[],2);
        end

        function y = getzd(obj,p)
            cp = [p(:,1)-obj.p0(1),p(:,2)-obj.p0(2)];
            index1 = abs(sqrt(sum(cp.^2,2))-obj.getRadius)<obj.geps;
            %             p = p(index,:);
            u12 = obj.getu12;
            u12 = [u12(2),-u12(1)];
            index2 = p*u12'-obj.p1*u12';
            if obj.direction > 0
                index2 = index2>=-obj.geps;
            else
                index2 = index2<=obj.geps;
            end
            y = bitand(index1,index2);
        end

        function y = getmirror(obj,varargin)
            y = obj;
            y.p1 = pmirror(y.p1,varargin{:});
            y.p2 = pmirror(y.p2,varargin{:});
            y.p0 = pmirror(y.p0,varargin{:});
            y.direction = -y.direction;
        end

        function y = getrotate(obj,varargin)
            y = obj;
            y.p1 = ext_protate2(y.p1,varargin{:});
            y.p2 = ext_protate2(y.p2,varargin{:});
            y.p0 = ext_protate2(y.p0,varargin{:});
        end

        function y = getshift(obj,varargin)
            y = obj;
            y.p1 = pshift(y.p1,varargin{:});
            y.p2 = pshift(y.p2,varargin{:});
            y.p0 = pshift(y.p0,varargin{:});
        end

        function y = getu1(obj)
            y = getVector(obj.p1-obj.p0)/getDistanceFromOrigin(obj.p1-obj.p0);
        end

        function y = getu12(obj)
            y = getVector(obj.p2-obj.p1)/getDistanceFromOrigin(obj.p2-obj.p1);
        end

        function y = getu2(obj)
            y = getVector(obj.p2-obj.p0)/getDistanceFromOrigin(obj.p2-obj.p0);
        end

        function y = getindent(obj,d)
            y = obj;
            u1 = y.getu1;
            u2 = y.getu2;
            y.p1 = y.p1 + d*u1;
            y.p2 = y.p2 + d*u2;
        end

        function y = getcta(obj,pname,side,radius)
            pname = rmspaces(pname);
            side = rmspaces(side);
            switch lower(pname)
                case 'p1'
                    u = obj.getu1;
                    switch lower(side)
                        case 'i'
                            y = obj.p1 - radius*u;
                        case 'o'
                            y = obj.p1 + radius*u;
                        otherwise
                            error('side muse be <<l>> or <<r>>.');
                    end
                case 'p2'
                    u = obj.getu2;
                    switch lower(side)
                        case 'i'
                            y = obj.p2 - radius*u;
                        case 'o'
                            y = obj.p2 + radius*u;
                        otherwise
                            error('side muse be <<i>> or <<o>>.');
                    end
                otherwise
                    error('pname must be <<p1>> or <<p2>>.');
            end
        end

        function y = getcpa(obj,pname,side,radius)
            pname = rmspaces(pname);
            side = rmspaces(side);
            switch lower(pname)
                case 'p1'
                    u = obj.getu1;
                    u = ext_protate2(u,pi/2);
                    switch lower(side)
                        case 'b'
                            y = obj.p1 - radius*u;
                        case 'f'
                            y = obj.p1 + radius*u;
                        otherwise
                            error('side muse be <<b>> or <<f>>.');
                    end
                case 'p2'
                    u = obj.getu2;
                    u = ext_protate2(u,pi/2);
                    switch lower(side)
                        case 'b'
                            y = obj.p2 - radius*u;
                        case 'f'
                            y = obj.p2 + radius*u;
                        otherwise
                            error('side muse be <<b>> or <<f>>.');
                    end
                otherwise
                    error('pname must be <<p1>> or <<p2>>.');
            end
        end

        function show(obj,label,scale)

            if nargin<3
                scale = 1;
            end

            p = [obj.p1;obj.getip;obj.p2];
            plot(p(:,1),p(:,2));
            n = size(p,1);
            switch rem(n,2)
                case 0
                    xp1 = p(n/2,:);
                    xp2 = p(n/2+1,:);
                    pc = (xp1 +  xp2)/2;
                case 1
                    xp1 = p(ceil(n/2)-1,:);
                    xp2 = p(ceil(n/2)+1,:);
                    pc = p(ceil(n/2),:);
            end
            vec = scale*(xp2-xp1)/norm(xp2-xp1);
            vec = ext_protate2(vec,pi/2);
            quiver(pc(1),pc(2),vec(1),vec(2),'color','r');
            text(pc(1)+vec(1)/2,pc(2)+vec(2)/2,label);
        end

        function y = getAngles(obj)
            % Angle of tangent lines at p1 and p2 w.r.t. x-axis, in [0, 2*pi)
            % y(1): tangent angle at p1
            % y(2): tangent angle at p2

            y = zeros(1,2);

            % radius angles from center to endpoints

            l_tmp = 1e-2/obj.getRadius;
            % tangent angles
            if obj.direction
                % CCW arc
                
                u = emdlab_g2d_rotatePoints(obj.p1.getVector, l_tmp,obj.p0.x,obj.p0.y);
                u = u - obj.p1.getVector;
%                 u = [-u(2),u(1)];
                y(1) = atan2(u(2),u(1));

                u = emdlab_g2d_rotatePoints(obj.p2.getVector, -l_tmp,obj.p0.x,obj.p0.y);
                u = u - obj.p2.getVector;
%                 u = [u(2),-u(1)];
                y(2) = atan2(u(2),u(1));

            else
                % CW arc

                u = emdlab_g2d_rotatePoints(obj.p1.getVector, -l_tmp,obj.p0.x,obj.p0.y);
                u = u - obj.p1.getVector;
%                 u = [u(2),-u(1)];
                y(1) = atan2(u(2),u(1));
                
                u = emdlab_g2d_rotatePoints(obj.p2.getVector, l_tmp,obj.p0.x,obj.p0.y);
                u = u - obj.p2.getVector;
%                 u = [-u(2),u(1)];
                y(2) = atan2(u(2),u(1));

            end

            y = mod(y, 2*pi);

        end


        function y = isID1(obj, pointID)
            y = obj.p1.id == pointID;
        end

        function y = isID2(obj, pointID)
            y = obj.p2.id == pointID;
        end

        function y = getPtr1(obj)
            y = obj.p1;
        end
        
        function y = getPtr2(obj)
            y = obj.p2;
        end

        function setPtr1(obj, pptr)
            obj.p1 = pptr;
        end

        function setPtr2(obj, pptr)
            obj.p2 = pptr;
        end

        function y = getTheta1Theta2(obj)

            tmp = obj.p1 - obj.p0;
            y(1) = tmp.getAngle;
            tmp = obj.p2 - obj.p0;
            y(2) = tmp.getAngle;
            y = rad2deg(y);

            if ~obj.direction
                y = fliplr(y);
            end

        end

        function y = isPointOnEdge(obj, p, tol)
            % Check whether point p lies on the arc interior within
            % tolerance tol, excluding the two endpoints.
            %
            % Inputs:
            %   p   : emdlab_g2d_point
            %   tol : distance tolerance
            %
            % Output:
            %   y   : logical true/false

            if nargin < 3
                tol = 1e-5;
            end

            R = obj.getRadius;

            % --- radius check: |dist(centre, p) - R| <= tol
            q = [p.x - obj.p0.x, p.y - obj.p0.y];
            r = norm(q);
            if abs(r - R) > tol || R <= tol
                y = false;   % not on the circle, or degenerate arc
                return;
            end

            % angular tolerance at this radius
            angTol = tol / R;

            % angle of the point relative to the start point, walked
            % in the direction of travel (CCW or CW)
            theta1 = atan2(obj.p1.y - obj.p0.y, obj.p1.x - obj.p0.x);
            thetaP = atan2(q(2), q(1));

            if obj.direction
                delta = mod(thetaP - theta1, 2*pi);      % CCW sweep
            else
                delta = mod(theta1 - thetaP, 2*pi);      % CW sweep
            end

            total = obj.getAngle;

            % --- exclude the endpoints
            if delta <= angTol || delta >= total - angTol
                y = false;
                return;
            end

            % --- point angularly inside the arc
            y = (delta <= total + angTol);
        end

        function y = getIndent(obj, shift)

            p0xy = obj.p0.getVector;
            u1 = obj.getu1;
            u2 = obj.getu2;

            p1xy = p0xy + (obj.getRadius - shift) * u1;
            p2xy = p0xy + (obj.getRadius - shift) * u2;

            y = emdlab_g2d_arc(emdlab_g2d_point(p0xy),emdlab_g2d_point(p1xy),emdlab_g2d_point(p2xy), obj.direction);

        end

        function xy = getPointProjection(obj, x, y)

            tmp = [x,y] - obj.p0.getVector;
            xy = obj.p0.getVector + obj.getRadius * tmp / norm(tmp);

        end

        function flg = isequal(obj, newObj)
            if isa(newObj, 'emdlab_g2d_arc')
                if (obj.p0.id == newObj.p0.id) && ((( obj.p1.id == newObj.p1.id) && (obj.p2.id == newObj.p2.id) && (obj.direction == newObj.direction)) || ...
                        (( obj.p1.id == newObj.p2.id) && (obj.p2.id == newObj.p1.id) && (obj.direction ~= newObj.direction)))
                    flg = true;
                    return;
                end                
            end
            flg = false;
        end

        function flg = hasIntersection(obj, newObj)
%             if isa(newObj, 'emdlab_g2d_segment')
%                 if ((obj.p0.id == newObj.p0.id) && (obj.p1.id == newObj.p1.id)) || ...
%                         ((obj.p0.id == newObj.p1.id) && (obj.p1.id == newObj.p0.id))
%                     flg = true;
%                     return;
%                 end    
%             else
%                 
%             end
            flg = true;
        end

    end

end
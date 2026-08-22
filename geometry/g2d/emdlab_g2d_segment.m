% EMDLAB: Electrical Machines Design Laboratory
% 2d segment class

classdef emdlab_g2d_segment < handle

    properties

        % initial and desination points
        p0 (1,1) emdlab_g2d_point;
        p1 (1,1) emdlab_g2d_point;

        % properties needed for wireframe mesh
        Nnodes (1,1) double {mustBeInteger, mustBeNonnegative} = 2;
        isSetNnodes (1,1) logical = true;
        maxLength (1,1) double {mustBeNonnegative} = 1;
        isSetMaxLength (1,1) logical = false;
        L1 (1,1) double {mustBeNonnegative} = 1;
        L2 (1,1) double {mustBeNonnegative} = 1;
        isSetL1L2 (1,1) logical = false;

    end

    methods
        %% constructor and destructor
        function obj = emdlab_g2d_segment(p0, p1)
            obj.p0 = p0;
            obj.p1 = p1;
        end

        function y = getLength(obj)
            y = getRadialLength(obj.p1 - obj.p0);
        end

        function y = getUnitVector(obj)
            y = getVector(obj.p1 - obj.p0);
            y = y/norm(y);
        end

        %% wireframe mesh functions
        function  nodes = getMeshNodesMinimal(obj)

%             nodes = zeros(obj.Nnodes,2);
%             nodes(:,1) = [obj.p0.x; obj.p1.x];
%             nodes(:,2) = [obj.p0.y; obj.p1.y];
            nodes = obj.getMeshNodes;

        end

        function nodes = getMeshNodes(obj)

            if obj.isSetNnodes

                nodes = zeros(obj.Nnodes,2);
                nodes(:,1) = linspace(obj.p0.x, obj.p1.x, obj.Nnodes);
                nodes(:,2) = linspace(obj.p0.y, obj.p1.y, obj.Nnodes);

            elseif obj.isSetMaxLength

                Nn = max(round(getRadialLength(obj.p1 - obj.p0) / obj.maxLength), 3);
                nodes = zeros(Nn,2);
                nodes(:,1) = linspace(obj.p0.x, obj.p1.x, Nn);
                nodes(:,2) = linspace(obj.p0.y, obj.p1.y, Nn);

            else

                % suppose we have n cords, n+1 points
                % [L1, a*L1, a^2*L1, ...., a^(n-1)*L1], a^(n-1)*L1 = L2
                L = obj.getLength;
                n = 2;
                a = nthroot(obj.L2/obj.L1,n-1);
                errOld = abs(L/obj.L1 - sum(a.^(0:n-1)));
                while true
                    n = n + 1;
                    a = nthroot(obj.L2/obj.L1,n-1);
                    errNew = abs(L/obj.L1 - sum(a.^(0:n-1)));
                    if errNew > errOld
                        n = n - 1;
                        break;
                    end
                    errOld = errNew;
                end

                err_fcn = @(x) sum(x.^(0:n-1))-L/obj.L1;
                a = fzero(err_fcn,1);
                u = obj.getUnitVector();
                nodes = zeros(n+1,2);
                nodes(:,1) = [obj.p0.x, obj.p0.x + cumsum(obj.L1 * u(1)* a.^(0:n-2)), obj.p1.x];
                nodes(:,2) = [obj.p0.y, obj.p0.y + cumsum(obj.L1 * u(2)* a.^(0:n-2)), obj.p1.y];

            end

        end

        function setNnodes(obj, n)

            obj.isSetMaxLength = false;
            obj.isSetNnodes = true;
            obj.isSetL1L2 = false;

            obj.Nnodes = n;

        end

        function setMaxLength(obj, mL)

            obj.isSetMaxLength = true;
            obj.isSetNnodes = false;
            obj.isSetL1L2 = false;

            obj.maxLength = mL;
            obj.p0.meshSize = mL;
            obj.p1.meshSize = mL;

        end

        function setL1L2(obj, newL1, newL2)

            if (newL1 + newL2) > obj.getLength()
                obj.setNnodes(2);
                return;
            end

            obj.isSetMaxLength = false;
            obj.isSetNnodes = false;
            obj.isSetL1L2 = true;

            obj.L1 = newL1;
            obj.L2 = newL2;

        end

        function y = getCenter(obj)

            y = [obj.p0.x + obj.p1.x, obj.p0.y + obj.p1.y]/2;

        end

        function y = getAngles(obj)
            % y(1): Angle of vector P0 -> P1
            % y(2): Angle of vector P1 -> P0

            y = zeros(1, 2);

            % Vector 1: P0 to P1
            dy1 = obj.p1.y - obj.p0.y;
            dx1 = obj.p1.x - obj.p0.x;
            y(1) = atan2(dy1, dx1);

            % Vector 2: P1 to P0 (This is mathematically just the angle of vector 1 + pi)
            % You don't need to recalculate atan2 here if you are confident in y(1)
            y(2) = atan2(-dy1, -dx1);

            y = mod(y, 2*pi);

        end

        function y = isID1(obj, pointID)
            y = obj.p0.id == pointID;
        end

        function y = isID2(obj, pointID)
            y = obj.p1.id == pointID;
        end

        function y = getPtr1(obj)
            y = obj.p0;
        end

        function y = getPtr2(obj)
            y = obj.p1;
        end

        function setPtr1(obj, pptr)
            obj.p0 = pptr;
        end

        function setPtr2(obj, pptr)
            obj.p1 = pptr;
        end

        function y = isPointOnEdge(obj, p, tol)
            % Check whether point p lies on the segment interior
            % within tolerance tol, excluding endpoints.
            %
            % Inputs:
            %   p   : emdlab_g2d_point
            %   tol : distance tolerance
            %
            % Output:
            %   y   : logical true/false

            if nargin < 3
                tol = 1e-9;
            end

            a = [obj.p0.x, obj.p0.y];
            b = [obj.p1.x, obj.p1.y];
            q = [p.x, p.y];

            ab = b - a;
            aq = q - a;

            L2_ = dot(ab, ab);
            if L2_ <= tol^2
                % Degenerate segment
                y = false;
                return;
            end

            % Projection parameter on the infinite line
            t = dot(aq, ab) / L2_;

            % Exclude endpoints with tolerance measured along the segment
            L = sqrt(L2_);
            if t < tol / L || t > 1 - tol / L
                y = false;
                return;
            end

            % Perpendicular distance from point to the line
            qproj = a + t * ab;
            d = norm(q - qproj);

            y = (d <= tol);
        end

        function y = getIndent(obj, shift)

            u = obj.getUnitVector;
            v = [-u(2),u(1)];
            xy0 = obj.p0.getVector;
            xy1 = obj.p1.getVector;

            xy0 = xy0 + shift * v;
            xy1 = xy1 + shift * v;

            p0ptr = emdlab_g2d_point(xy0(1),xy0(2));
            p1ptr = emdlab_g2d_point(xy1(1),xy1(2));
            y = emdlab_g2d_segment(p0ptr,p1ptr);

        end

        function xy = getPointProjection(obj, x, y)

            tmp = [x,y] - obj.p0.getVector;
            u = obj.getUnitVector;
            xy = obj.p0.getVector + (u*tmp') * u;

        end

        function flg = isequal(obj, newObj)
            if isa(newObj, 'emdlab_g2d_segment')
                if ((obj.p0.id == newObj.p0.id) && (obj.p1.id == newObj.p1.id)) || ...
                        ((obj.p0.id == newObj.p1.id) && (obj.p1.id == newObj.p0.id))
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
%             elseif isa(newObj, 'emdlab_g2d_arc')
%                 
%             end
            flg = true;
        end

    end

end
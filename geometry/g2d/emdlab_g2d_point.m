% EMDLAB: Electrical Machines Design Laboratory
% 2d point class

classdef emdlab_g2d_point  < handle & matlab.mixin.Copyable & emdlab_g2d_constants

    properties

        % x and y coordinates of the point
        x (1,1) double = 0;
        y (1,1) double = 0;

        % the object tag
        tag (1,:) char = '';

        % tag of edges that used this point for their constrcution
        tags (1,:) string;

        % mesh size: mesh size at this point, imporant in Gmsh
        meshSize (1,1) double = 1;

    end

    methods
        %% constructor: you can define the point with (x,y) input or ([x,y])
        function obj = emdlab_g2d_point(varargin)

            if nargin == 1

                obj.x = varargin{1}(1);
                obj.y = varargin{1}(2);

            elseif nargin == 2

                obj.x = varargin{1};
                obj.y = varargin{2};

            end

        end

        %% set and get coordinates
        function setX(obj, newX)
            obj.x = newX;
        end

        function setY(obj, newY)
            obj.y = newY;
        end

        function setCoordinates(obj, newX, newY)
            obj.x = newX;
            obj.y = newY;
        end

        function val = getX(obj)
            val = obj.x;
        end

        function val = getY(obj)
            val = obj.y;
        end

        % get vector format of the point
        function vec = getVector(obj)

            vec = [obj.x, obj.y];

        end

        % set in vector format ([x,y])
        function setVector(obj, vec)

            obj.x = vec(1);
            obj.y = vec(2);

        end

        function u = getUnitVector(obj)

            u = obj.getVector/obj.getRadialLength;

        end

        function normalize(obj)

            dfo = obj.getDistanceFromOrigin;
            obj.x = obj.x / dfo;
            obj.y = obj.y / dfo;

        end

        %% operation functions
        % add two points
        function p2 = plus(obj, p1)

            p2 = emdlab_g2d_point(obj.x + p1.x, obj.y + p1.y);

        end

        % subtract two point
        function p2 = minus(obj, p1)

            p2 = emdlab_g2d_point(obj.x - p1.x, obj.y - p1.y);

        end

        function p2 = mtimes(lInput, rInput)

            if isa(lInput, 'emdlab_g2d_point') && isa(rInput, 'emdlab_g2d_point')

                p2 = emdlab_g2d_point(0,0);
                p2.x = lInput.x * rInput.x - lInput.y * rInput.y;
                p2.y = lInput.x * rInput.y - lInput.y * rInput.x;

            elseif isa(lInput, 'emdlab_g2d_point') && isscalar(rInput)

                p2 = emdlab_g2d_point(rInput*lInput.x, rInput*lInput.y);

            elseif isa(rInput, 'emdlab_g2d_point') && isscalar(lInput)

                p2 = emdlab_g2d_point(lInput*rInput.x, lInput*rInput.y);

            else
                error('Unsupported multiplication type.');
            end

        end

        %% distance functions
        % distance of the point from origin
        function d = getRadialLength(obj)
            d = sqrt(obj.x.^2 + obj.y.^2);
        end

        function d = norm(obj)
            d = obj.getRadialLength;
        end

        function d = getDistanceFromOrigin(obj)
            d = sqrt(obj.x.^2 + obj.y.^2);
        end

        function d = getDistanceFromPoint(obj, p)
            d = sqrt((obj.x - p.x).^2 + (obj.y - p.y).^2);
        end        

        function d = getDistanceFromLine(obj, l)

            arguments
                obj (1,1) emdlab_g2d_point;
                l (1,1) emdlab_g2d_line;
            end

            p0p = obj - l.p0;
            n = p0p - p0p.dot(l.u) * l.u;
            d = n.getDistanceFromOrigin;

        end

        function setDistanceFromPoint(obj, p, d)

            u12 = obj - p;
            u12.normalize;
            obj.x = p.x + d*u12.x;
            obj.y = p.y + d*u12.y;

        end

        function setDistanceFromLine(obj, l, d)

            u12 = obj - l.p0;
            vp = l.u.dot(u12) * l.u;
            vn = u12 - vp;
            vn.normalize;
            newPoint = l.p0 + vp + d*vn;
            obj.setCoordinates(newPoint.x, newPoint.y);

        end

        %% inner and outer products
        % outer product of two points
        function val = outerProduct(obj, p)
            val = obj.x * p.y - obj.y * p.x;
        end

        function val = cross(obj, p1)
            val = obj.x * p1.y - obj.y * p1.x;
        end

        % inner product of two points
        function val = innerProduct(obj, p)
            val = obj.x * p.x + obj.y * p.y;
        end

        function val = dot(obj, p1)
            val = obj.x * p1.x + obj.y * p1.y;
        end

        %% rotate functions
        function rotateAroundOrigin(obj, angle)

            obj.setVector(ext_protate2(obj.getVector, angle));

        end

        function rotateAroundPoint(obj, p0, angle)

            obj.setVector(ext_protate2(obj.getVector, angle, p0.getVector));

        end

        function p = getRotateAroundOrigin(obj, angle)

            p = emdlab_g2d_point(0,0);
            p.setVector(ext_protate2(obj.getVector, angle));

        end

        function p = getRotateAroundPoint(obj, p0, angle)

            p = emdlab_g2d_point(0,0);
            p.setVector(ext_protate2(obj.getVector, angle, p0.getVector));

        end

        %% shift functions
        function shift(obj, xsh, ysh)

            obj.x = obj.x + xsh;
            obj.y = obj.y + ysh;

        end

        function shiftX(obj, xsh)

            obj.x = obj.x + xsh;

        end

        function shiftY(obj, ysh)

            obj.y = obj.y + ysh;

        end

        function p = getShift(obj, varargin)

            p = copy(obj);
            p.shift(varargin{:});

        end

        function p = getShiftX(obj, varargin)

            p = copy(obj);
            p.shiftX(varargin{:});

        end

        function p = getShiftY(obj, varargin)

            p = copy(obj);
            p.shiftY(varargin{:});

        end

        %% move functions
        function moveX(obj, deltaX)
            obj.x = obj.x + deltaX;
        end

        function moveY(obj, deltaY)
            obj.y = obj.y + deltaY;
        end

        function moveXY(obj, deltaX, deltaY)
            obj.x = obj.x + deltaX;
            obj.y = obj.y + deltaY;
        end

        function moveRadial(obj, deltaR)
            u = obj.getVector / obj.getRadialLength;
            obj.x = obj.x + deltaR * u(1);
            obj.y = obj.y + deltaR * u(2);
        end        

        function newObj = getMoveX(obj, varargin)
            newObj = copy(obj);
            newObj.moveX(varargin{:});
        end

        function newObj = getMoveY(obj, varargin)
            newObj = copy(obj);
            newObj.moveY(varargin{:});
        end

        function newObj = getMoveXY(obj, varargin)
            newObj = copy(obj);
            newObj.moveXY(varargin{:});
        end

        function newObj = getMoveRadial(obj, varargin)
            newObj = copy(obj);
            newObj.moveRadial(varargin{:});
        end

        function y = getAngle(obj)
            y = mod(atan2(obj.y, obj.x), 2*pi);
        end

    end

end

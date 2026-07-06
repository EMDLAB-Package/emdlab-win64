classdef emdlab_m3d_cs 

    properties
        origin = [0,0,0];
        ux = [1,0,0];
        uy = [0,1,0];
        uz = [0,0,1];
    end

    methods

        function obj = emdlab_m3d_cs
        end

        function obj = setCoordinateSystem(obj, xAxis, yPoint, originPoint)

            if nargin<3
                originPoint = [0,0,0];
            end

            obj.origin = originPoint;

            obj.ux = xAxis - originPoint;
            obj.ux = obj.ux /norm(obj.ux);

            tmp = yPoint - originPoint;

            obj.uz = cross(obj.ux, tmp);
            obj.uz = obj.uz /norm(obj.uz);

            obj.uy = cross(obj.uz, obj.ux);

        end

        function xyz = getCoordinateSystem(obj)
            
        end

    end

end
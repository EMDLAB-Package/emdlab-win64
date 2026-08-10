% EMDLAB: Electrical Machines Design Laboratory
% 2d edge class
% an edge can be:
% a segment
% an arc
% a spline

classdef emdlab_g2d_edge  < handle & matlab.mixin.Copyable & emdlab_g2d_constants

    properties

        % pointer to edge object type
        ptr (1,1);

        % the object tag
        tag (1,:) char = '';

        % tag of loops that used this edge for their constrcution
        tags (1,:) string;

    end

    methods

        function flag = isSegment(obj)
            if isa(obj.ptr, 'emdlab_g2d_segment')
                flag = true;
            else
                flag = false;
            end
        end

        function flag = isArc(obj)
            if isa(obj.ptr, 'emdlab_g2d_arc')
                flag = true;
            else
                flag = false;
            end
        end

        function flag = isSpline(obj)
            if isa(obj.ptr, 'emdlab_g2d_spline')
                flag = true;
            else
                flag = false;
            end
        end

    end

end
% EMDLAB: Electrical Machines Design Laboratory
% emdlab data base class for 2d geometries

classdef emdlab_g2d_db < handle

    properties

        % points
        points (1,:) emdlab_g2d_point;

        % edges
        edges (1,:) emdlab_g2d_edge;

        % loops
        loops (1,:) emdlab_g2d_loop;

        % faces
        faces (1,:) emdlab_g2d_face;

    end

    properties (SetAccess=protected)
        % python path
        pyPath = "";
    end

    methods
        %% constructor and destructor
        function obj = emdlab_g2d_db()

            % set python path
            p = pyenv;
            if p.Executable ~= ""
                obj.setPyPath(p.Executable);
            end

        end

        function setPyPath(obj, filePath)
            if ~isfile(filePath)
                error('There is no file in specified path.');
            end
            obj.pyPath = string(filePath);
            obj.pyPath = replace(obj.pyPath, '\', '\\');
        end

        %% point methods
        % adding a new point to data base
        % this function returns point index and point handle
        function varargout = addPoint(obj, varargin)

            % check the varargin type
            if numel(varargin) == 2
                x = varargin{1};
                y = varargin{2};
            elseif numel(varargin) == 1
                if isa(varargin{1},'emdlab_g2d_point')
                    x = varargin{1}.x;
                    y = varargin{1}.y;
                elseif isvector(varargin{1}) && isnumeric(varargin{1}) && (length(varargin{1}) == 2)
                    x = varargin{1}(1);
                    y = varargin{1}(2);
                else
                    throw(MException('', 'Wrong input type, it must be <emdlab_g2d_point> type.'));
                end
            else
                throw(MException('', 'Wrong number of input arguments.'));
            end

            % check for existance of already defined point in the same location
            for i = 1:numel(obj.points)

                if norm(obj.points(i).getVector() - [x,y]) < 1e-6

                    pointHandle = obj.points(i);
                    if nargout == 1
                        varargout{1} = i;
                    elseif nargout == 2
                        varargout{1} = i;
                        varargout{2} = pointHandle;
                    elseif nargout > 2
                        error('The number of output arguments is too high.');
                    end
                    return;

                end

            end

            % get an instance of point class
            pointHandle = emdlab_g2d_point(x,y);
            obj.points(end+1) = pointHandle;

            % generate point tag
            pointIndex = numel(obj.points);
            pointHandle.tag = ['p', num2str(pointIndex)];

            if nargout == 1
                varargout{1} = pointIndex;
            elseif nargout == 2
                varargout{1} = pointIndex;
                varargout{2} = pointHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function removePoint(obj, pIndex)

            % remove connected edges to this point firt
            for i = 1:numel(pIndex.tags)
            end
            obj.points(pIndex) = [];
        end

        function pointHandle = getPointHandleByTag(obj, pTag)

            % check for existance of already defined point in data base
            for i = 1:numel(obj.points)

                if strcmp(obj.points(i).tag,pTag)

                    pointHandle = obj.points(i);
                    return;

                end

            end

            error('Point was not found.');

        end

        function pointIndex = getPointIndexByTag(obj, pTag)

            % check for existance of already defined point in data base
            for i = 1:numel(obj.points)

                if strcmp(obj.points(i).tag,pTag)

                    pointIndex = i;
                    return;

                end

            end

            error('Point was not found.');

        end

        function pointIndex = getPointIndexByCoordinates(obj, x, y)

            if ~numel(obj.points)
                error('There is no defined point.');
            end

            % this function returns index of closed point to x and y coordinates
            minDistance = inf;
            for i = 1:numel(obj.points)

                distance = norm(obj.points(i).getVector - [x,y]);
                if distance < minDistance
                    pointIndex = i;
                    minDistance = distance;
                end

            end

        end

        function [x,y] = getPointCoordinates(obj, pIndex)
            x = obj.points(pIndex).x;
            y = obj.points(pIndex).y;
        end

        function alignPointsAlongYAxis(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).x = obj.points(varargin{1}).x;
            end

        end

        function alignPointsAlongXAxis(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsAlongRAxis(obj, varargin)

            u_ref = obj.points(varargin{1}).getUnitVector;
            for i = 2:numel(varargin)
                r_i = obj.points(varargin{i}).getDistanceFromOrigin;
                obj.points(varargin{i}).setCoordinates(r_i*u_ref(1), r_i*u_ref(2));
            end

        end

        function alignPointsAlongTAxis(obj, varargin)

            r_ref = obj.points(varargin{1}).getDistanceFromOrigin;
            for i = 2:numel(varargin)
                u_i = obj.points(varargin{i}).getUnitVector;
                obj.points(varargin{i}).setCoordinates(r_ref*u_i(1), r_ref*u_i(2));
            end

        end

        % align points functions
        function alignPointsXX(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).x = obj.points(varargin{1}).x;
            end

        end

        function alignPointsXY(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsXR(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsXT(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsYX(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsYY(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsYR(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsYT(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsRX(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsRY(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsRR(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsRT(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsTX(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsTY(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsTR(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsTT(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function setPointDistanceFromPoint(obj, pIndex, x, y, distance)

            if abs(norm([obj.points(pIndex).x - x, obj.points(pIndex).y - y])) < 1e-6
                error(' Points are on each other.');
            end

            u = obj.points(pIndex) - emdlab_g2d_point(x,y);
            u.normalize();
            obj.points(pIndex).x = x + distance * u.x;
            obj.points(pIndex).y = y + distance * u.y;

        end

        function setPointDistanceFromP0ULine(obj, pIndex, x, y, ux, uy, distance)
            obj.points(pIndex).setDistanceFromLine(emdlab_g2d_line(x,y,ux,uy),distance);
        end

        function str = getPointsXCoordinatesForMaxwell(obj, pIndex, unit)

            if nargin<3, unit = 'mm'; end
            str = "[";
            for pi = pIndex
                str = str + sprintf('%.12g', obj.points(pi).x) + ",";
            end
            str = char(str);
            str(end) = "]";
            str = str + " " + unit;

        end

        function str = getPointsYCoordinatesForMaxwell(obj, pIndex, unit)

            if nargin<3, unit = 'mm'; end
            str = "[";
            for pi = pIndex
                str = str + sprintf('%.12g', obj.points(pi).y) + ",";
            end
            str = char(str);
            str(end) = "]";
            str = str + " " + unit;

        end

        function str = getEdgesAnglesForMaxwell(obj, eIndex)

            str = "[";
            for ei = eIndex
                if isa(obj.edges(ei).ptr, 'emdlab_g2d_arc')
                    if obj.edges(ei).ptr.direction
                        str = str + sprintf('%.12g', real(obj.edges(ei).ptr.getAngleDegree)) + ",";
                    else
                        str = str + sprintf('-%.12g', real(obj.edges(ei).ptr.getAngleDegree)) + ",";
                    end
                else
                    str = str + "0,";
                end
            end
            str = char(str);
            str(end) = "]";
            str = str + " deg";

        end

        %% edge methods
        % adding a new segment to data base
        % this function returns edge index and edge handle
        function varargout = addSegment(obj, p0Index, p1Index)

            % check for char inputs
            if ischar(p0Index)
                p0Index = obj.getPointIndexByTag(p0Index);
            end

            if ischar(p1Index)
                p1Index = obj.getPointIndexByTag(p1Index);
            end

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;
            % set pointer class to segment
            edgeHandle.ptr = emdlab_g2d_segment(obj.points(p0Index),obj.points(p1Index));
            obj.edges(end+1) = edgeHandle;

            % generate edge tag
            edgeIndex = numel(obj.edges);
            edgeHandle.tag = ['e', num2str(edgeIndex)];
            edgeHandle.ptr.tag = edgeHandle.tag;

            % add edge tag to connected points
            obj.points(p0Index).tags(end+1) = edgeHandle.tag;
            obj.points(p1Index).tags(end+1) = edgeHandle.tag;

            if nargout == 1
                varargout{1} = edgeIndex;
            elseif nargout == 2
                varargout{1} = edgeIndex;
                varargout{2} = edgeHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % adding a new spline to data base
        % this function returns edge index and edge handle
        function varargout = addSpline(obj, ptsIndex)

            % number of points
            Np = numel(ptsIndex);
            pts = repmat(emdlab_g2d_point,1,Np);

            if isa(ptsIndex, 'cell')

                for i = 1:Np
                    if ischar(ptsIndex{i})
                        pts(i) = obj.points(obj.getPointIndexByTag(ptsIndex{i}));
                    else
                        pts(i) = obj.points(ptsIndex{i});
                    end
                end

            else

                for i = 1:Np
                    pts(i) = obj.points(ptsIndex(i));
                end

            end

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;

            % set pointer class to segment
            edgeHandle.ptr = emdlab_g2d_spline(pts);
            obj.edges(end+1) = edgeHandle;

            % generate edge tag
            edgeIndex = numel(obj.edges);
            edgeHandle.tag = ['e', num2str(edgeIndex)];
            edgeHandle.ptr.tag = edgeHandle.tag;

            if nargout == 1
                varargout{1} = edgeIndex;
            elseif nargout == 2
                varargout{1} = edgeIndex;
                varargout{2} = edgeHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % add a new segment by direct coordinates passing
        function varargout = addSegmentByCoordinates(obj, x1, y1, x2, y2)

            p1Index = obj.addPoint(x1,y1);
            p2Index = obj.addPoint(x2,y2);

            if nargout == 0
                obj.addSegment(p1Index,p2Index);
            elseif nargout == 1
                varargout{1} = obj.addSegment(p1Index,p2Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegment(p1Index,p2Index);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % add a new segment by direct coordinates passing
        function varargout = addSplineByCoordinates(obj, x, y)

            % check feasibility
            if length(x) ~= length(y)
                error('x and y vectors must have the same length.');
            end

            Nx = length(x);
            ptsIndex = zeros(1,Nx);
            for i = 1:Nx
                ptsIndex(i) = obj.addPoint(x(i),y(i));
            end

            if nargout == 0
                obj.addSpline(ptsIndex);
            elseif nargout == 1
                varargout{1} = obj.addSpline(ptsIndex);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSpline(ptsIndex);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % adding a new arc to data base
        function varargout = addArc(obj, p0Index, p1Index, p2Index, direction)

            % check for char inputs
            if ischar(p0Index)
                p0Index = obj.getPointIndexByTag(p0Index);
            end

            if ischar(p1Index)
                p1Index = obj.getPointIndexByTag(p1Index);
            end

            if ischar(p2Index)
                p2Index = obj.getPointIndexByTag(p2Index);
            end

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;
            % set pointer class to arc
            edgeHandle.ptr = emdlab_g2d_arc(obj.points(p0Index),obj.points(p1Index),obj.points(p2Index), direction);
            obj.edges(end+1) = edgeHandle;

            % generate edge tag
            edgeIndex = numel(obj.edges);
            edgeHandle.tag = ['e', num2str(edgeIndex)];
            edgeHandle.ptr.tag = edgeHandle.tag;

            % add edge to point tags
            obj.points(p0Index).tags(end+1) = edgeHandle.tag;
            obj.points(p1Index).tags(end+1) = edgeHandle.tag;
            obj.points(p2Index).tags(end+1) = edgeHandle.tag;

            if nargout == 1
                varargout{1} = edgeIndex;
            elseif nargout == 2
                varargout{1} = edgeIndex;
                varargout{2} = edgeHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addArcByCoordinates(obj, x1, y1, x2, y2, x3, y3, direction)

            p1Index = obj.addPoint(x1,y1);
            p2Index = obj.addPoint(x2,y2);
            p3Index = obj.addPoint(x3,y3);

            if nargout == 0
                obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout == 1
                varargout{1} = obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addArcByCoordinatesCPA(obj, x1, y1, x2, y2, arcAngle)

            p1Index = obj.addPoint(x1,y1);
            p2Index = obj.addPoint(x2,y2);
            [x3,y3] = emdlab_g2d_rotatePointsXY(x2,y2,arcAngle,x1,y1);
            p3Index = obj.addPoint(x3,y3);
            direction = arcAngle > 0;

            if nargout == 0
                obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout == 1
                varargout{1} = obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % get edge handle
        function edgeHandle = getEdgeHandleByIndex(obj, eIndex)

            % check for existance of already defined edge in data base
            if eIndex <= numel(obj.edges)
                edgeHandle = obj.edges(eIndex).ptr;
            else
                error('Edge was not found.');
            end

        end

        function edgeHandle = getEdgeHandleByTag(obj, eTag)

            % check for existance of already defined edge in data base
            for i = 1:numel(obj.edges)

                if strcmp(obj.edges(i).tag,eTag)

                    edgeHandle = obj.edges(i).ptr;
                    return;

                end

            end

            error('Edge was not found.');

        end

        function edgeIndex = getEdgeIndexByTag(obj, eTag)

            % check for existance of already defined edge in data base
            for i = 1:numel(obj.edges)

                if strcmp(obj.edges(i).tag,eTag)

                    edgeIndex = i;
                    return;

                end

            end

            error('Edge was not found.');

        end

        function pIndex = getEdgeStartPointIndex(obj, eIndex)

            % pointer to edge
            eptr = obj.edges(eIndex).ptr;
            switch class(eptr)
                case 'emdlab_g2d_segment'
                    pIndex = obj.getPointIndexByTag(eptr.p0.tag);
                case 'emdlab_g2d_arc'
                    pIndex = obj.getPointIndexByTag(eptr.p1.tag);
                case 'emdlab_g2d_spline'
                    pIndex = obj.getPointIndexByTag(eptr.pts(1).tag);
            end

        end

        function pIndex = getEdgeEndPointIndex(obj, eIndex)

            % pointer to edge
            eptr = obj.edges(eIndex).ptr;
            switch class(eptr)
                case 'emdlab_g2d_segment'
                    pIndex = obj.getPointIndexByTag(eptr.p1.tag);
                case 'emdlab_g2d_arc'
                    pIndex = obj.getPointIndexByTag(eptr.p2.tag);
                case 'emdlab_g2d_spline'
                    pIndex = obj.getPointIndexByTag(eptr.pts(end).tag);
            end

        end

        function splitSegment(obj, eIndex)

            edgeHandle = obj.edges(eIndex).ptr;
            tmp = edgeHandle.getCenter;
            p = obj.addPoint(tmp(1),tmp(2));
            p2 = edgeHandle.p1;
            edgeHandle.p1 = obj.points(p);
            obj.addSegment(p,obj.addPoint(p2));

        end

        function newEdgeIndex = splitArc(obj, eIndex)

            edgeHandle = obj.edges(eIndex).ptr;
            tmp = edgeHandle.p0.getVector;
            tmp = emdlab_g2d_rotatePoints(edgeHandle.p1.getVector,edgeHandle.getSignedAngle/2, tmp(1),tmp(2));
            p = obj.addPoint(tmp(1),tmp(2));
            p2 = edgeHandle.p2;
            edgeHandle.p2 = obj.points(p);
            newEdgeIndex = obj.addArc(obj.addPoint(edgeHandle.p0.getVector),p,obj.addPoint(p2),edgeHandle.direction);

        end

        % edge extensions to draw complex geometries
        function varargout = extendSegmentBySegmentUpToPoint(obj, eIndex, x, y, seIndex)

            % set default start/end index as end
            if nargin < 5, seIndex = 1; end

            % get segment edge handle
            edgeHandle = obj.edges(eIndex).ptr;

            % check start or index: 0 when we exten inward, 0: outward extension
            if seIndex == 0
                v2 = edgeHandle.p0.getVector;
                v1 = [x,y];
            else
                v1 = edgeHandle.p1.getVector;
                v2 = [x,y];
            end

            if nargout == 0
                obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 1
                varargout{1} = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendSegmentBySegment(obj, eIndex, extAngle, extAmplitude, seIndex)

            % set default start/end index as end
            if nargin < 5, seIndex = 1; end

            % get segment edge handle
            edgeHandle = obj.edges(eIndex).ptr;
            u = edgeHandle.getUnitVector;
            u = emdlab_g2d_rotatePoints(u,extAngle);

            % check start or index: 0 when we exten inward, 0: outward extension
            if seIndex == 0
                v2 = edgeHandle.p0.getVector;
                v1 = v2 + extAmplitude * u;
            else
                v1 = edgeHandle.p1.getVector;
                v2 = v1 + extAmplitude * u;
            end

            if nargout == 0
                obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 1
                varargout{1} = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendSegmentByTangentArc(obj, eIndex, extRadius, extAngle, seIndex)

            % set default start/end index as end
            if nargin < 5, seIndex = 1; end

            edgeHandle = obj.edges(eIndex).ptr;
            u = edgeHandle.getUnitVector;
            u = emdlab_g2d_rotatePoints(u,pi/2);

            % check start or index: 0 when we exten inward, 0: outward extension
            if seIndex == 0
                v2 = edgeHandle.p0.getVector;
                c = v2 + extRadius * u;
                v1 = emdlab_g2d_rotatePoints(v2, -sign(extRadius)*extAngle, c(1), c(2));
            else
                v1 = edgeHandle.p1.getVector;
                c = v1 + extRadius * u;
                v2 = emdlab_g2d_rotatePoints(v1, sign(extRadius)*extAngle, c(1), c(2));
            end

            if nargout == 0
                obj.addArcByCoordinates(c(1),c(2),v1(1),v1(2),v2(1),v2(2),extRadius>0);
            elseif nargout == 1
                varargout{1} = obj.addArcByCoordinates(c(1),c(2),v1(1),v1(2),v2(1),v2(2),extRadius>0);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArcByCoordinates(c(1),c(2),v1(1),v1(2),v2(1),v2(2),extRadius>0);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendSegmentByArc(obj, eIndex, xc, yc, extAngle, seIndex)

            % set default start/end index as end
            if nargin < 6, seIndex = 1; end

            % get edge handle
            edgeHandle = obj.edges(eIndex).ptr;

            % check start or end point of the segment for extension
            if seIndex == 0
                v1 = edgeHandle.p0.getVector;
            else
                v1 = edgeHandle.p1.getVector;
            end

            v2 = emdlab_g2d_rotatePoints(v1, extAngle, xc, yc);
            if nargout == 0
                obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout == 1
                varargout{1} = obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendArcBySegment(obj, eIndex, extAngle, extAmplitude, seIndex)

            % set default start/end index as end
            if nargin < 5, seIndex = 1; end

            % get arc edge handle
            edgeHandle = obj.edges(eIndex).ptr;

            % check start/end point index: 0 when we extend inward, 1: for outward extension
            if seIndex == 0

                u = edgeHandle.getu1;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,pi/2+extAngle);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 - extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,-pi/2+extAngle);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 - extAmplitude * u;
                end

            elseif seIndex == 1

                u = edgeHandle.getu2;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,pi/2+extAngle);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,-pi/2+extAngle);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                end

            else

                error('Start/end index must be 0 or 1.');

            end

            if nargout == 0
                obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 1
                varargout{1} = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendArcByTangentSegment(obj, eIndex, seIndex, extAmplitude)

            edgeHandle = obj.edges(eIndex).ptr;

            % check start or index: 0 when we exten inward, 0: outward extension
            if seIndex == 0

                u = edgeHandle.getu1;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,-pi/2);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,pi/2);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 + extAmplitude * u;
                end

            else

                u = edgeHandle.getu2;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,pi/2);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,-pi/2);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                end

            end

            if nargout == 0
                obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 1
                varargout{1} = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendArcByPerpendicularSegment(obj, eIndex, seIndex, extAmplitude)

            % get arc handle
            edgeHandle = obj.edges(eIndex).ptr;

            % check start or index: 0 when we exten inward, 0: outward extension
            if seIndex == 0

                u = -edgeHandle.getu1;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,-pi/2);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,pi/2);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 + extAmplitude * u;
                end

            else

                u = edgeHandle.getu2;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,pi/2);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,-pi/2);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                end

            end

            if nargout == 0
                obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 1
                varargout{1} = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendArcByArc(obj, eIndex, xc, yc, extAngle, seIndex)

            % set default start/end index as end
            if nargin < 6, seIndex = 1; end

            % get edge handle
            edgeHandle = obj.edges(eIndex).ptr;

            % check start or end point of the segment for extension
            if seIndex == 0
                v1 = edgeHandle.p1.getVector;
            else
                v1 = edgeHandle.p2.getVector;
            end

            v2 = emdlab_g2d_rotatePoints(v1, extAngle, xc, yc);
            if nargout == 0
                obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout == 1
                varargout{1} = obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function y = getEdgeLength(obj, eIndex)
            y = obj.edges(eIndex).ptr.getLength;
        end

        %% edge edits
        function indicesOfNewEdges = splitEdge(obj, edgeIndex, splitRatio)

            % splitRatio is a vector of ratios
            if any(splitRatio >= 1) || any(splitRatio <= 0)
                error('All split ratios must be between 0 and 1.');
            end

            if sum(splitRatio) > 1
                error('The summation of the split ratio mast be lower than one.');
            end

            % get edge pointer
            eptr = obj.edges(edgeIndex).ptr;

            % vector containing edge indicies
            indicesOfNewEdges = zeros(1,length(splitRatio)+1);
            indicesOfNewEdges(1) = edgeIndex;
            indicesOfNewPoints = zeros(1,length(splitRatio));

            % detect edge type and apply proper splits
            switch class(eptr)

                case 'emdlab_g2d_segment'

                    % get index of edge end point
                    p1Index = obj.getPointIndexByTag(eptr.p1.tag);

                    % calculate new points to split the edge
                    newp = eptr.p0.getVector;
                    vec = eptr.getUnitVector * eptr.getLength;
                    for i = 1:length(splitRatio)
                        newp = newp + splitRatio(i) * vec;
                        indicesOfNewPoints(i) = obj.addPoint(newp(1),newp(2));
                    end

                    % modify the first edge end point
                    eptr.p1 = obj.points(indicesOfNewPoints(1));

                    % add new edges
                    for i = 1:length(splitRatio)-1
                        indicesOfNewEdges(i+1) = obj.addSegment(indicesOfNewPoints(i),indicesOfNewPoints(i+1));
                    end

                    % add last edge
                    indicesOfNewEdges(end) = obj.addSegment(indicesOfNewPoints(end),p1Index);

            end

        end

        function removeEdge(obj, eIndex)

            % first remove all connected loops
            for lTag = obj.edges(eIndex).tags

            end

        end

        %% loop methods
        % adding a new loop to data base
        % this function returns loop index and loop handle
        function varargout = addLoop(obj, varargin)

            % get a loop class instance
            loopHandle = emdlab_g2d_loop;
            obj.loops(end+1) = loopHandle;

            % loop to assign edges and directions
            for i = 1:numel(varargin)

                for j = 1:numel(varargin{i})

                    if ischar(varargin{i}(j))

                        % get edge tag and remove spaces
                        edgeTag = strrep(varargin{i}(j), ' ', '');

                        if strcmpi(edgeTag(1), '-')
                            edgeIndex = obj.getEdgeIndexByTag(edgeTag(2:end));
                            loopHandle.addEdge(obj.edges(edgeIndex).ptr, false, -edgeIndex);
                        else
                            edgeIndex = obj.getEdgeIndexByTag(edgeTag);
                            loopHandle.addEdge(obj.edges(edgeIndex).ptr, true, edgeIndex);
                        end

                    else

                        edgeIndex = abs(varargin{i}(j));
                        loopHandle.addEdge(obj.edges(edgeIndex).ptr, varargin{i}(j)>0, varargin{i}(j));

                    end

                end

            end

            % generate loop tag
            loopIndex = numel(obj.loops);
            loopHandle.tag = ['l', num2str(loopIndex)];

            % add loop tag to connected edges
            for edgeIndex = abs(loopHandle.edgesIndexList)
                obj.edges(edgeIndex).tags(end+1) = loopHandle.tag;
            end

            if nargout == 1
                varargout{1} = loopIndex;
            elseif nargout == 2
                varargout{1} = loopIndex;
                varargout{2} = loopHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function removeLoop(obj, loopIndex)

            if ischar(loopIndex)
                loopIndex = obj.getLoopIndexByTag(loopIndex);
            end

            % first remove all connected faces to this loop
            for faceTag = obj.loops(loopIndex).tags
                obj.removeFace(faceTag);
            end

            obj.loops(loopIndex) = [];

        end

        function loopIndex = getLoopIndexByTag(obj, lTag)

            % check for existance of already defined loop in data base
            for i = 1:numel(obj.loops)

                if strcmp(obj.loops(i).tag,lTag)

                    loopIndex = i;
                    return;

                end

            end

            error('Loop was not found.');

        end

        %% face methods
        % adding a new face to data base
        % this function returns the face handle
        function varargout = addFace(obj, faceName, varargin)

            % get face class instance
            faceHandle = emdlab_g2d_face;
            faceHandle.tag = faceName;
            faceHandle.color = rand(1,3);
            obj.faces(end+1) = faceHandle;

            for i = 1:numel(varargin)

                for j = 1:numel(varargin{i})
                    faceHandle.addLoop(obj.loops(varargin{i}(j)));
                    obj.loops(varargin{i}(j)).tags(end+1) = faceName;
                end

            end

            if nargout == 1
                varargout{1} = faceHandle;
            elseif nargout > 1
                error('The number of output arguments is too hight');
            end

        end

        function removeFace(obj, faceIndex)

            if ischar(faceIndex)
                faceIndex = obj.getFaceIndexByTag(faceIndex);
            end
            obj.faces(faceIndex) = [];

        end

        function faceIndex = getFaceIndexByTag(obj, fTag)

            % check for existance of already defined face in data base
            for i = 1:numel(obj.faces)

                if strcmp(obj.faces(i).tag,fTag)

                    faceIndex = i;
                    return;

                end

            end

            error('Face was not found.');

        end

        function setFaceColor(obj, faceTag, R, G, B)
            obj.faces(obj.getFaceIndexByTag(faceTag)).color = [R,G,B]/255;
        end

        %% mesh generation methods
        % generate triangular mesh for geometry
        function m = generateMesh(obj, meshGenerator)

            % default mesh generator
            if nargin<2
                meshGenerator = 'mm';
            end

            if strcmpi(meshGenerator, 'gmsh')
                obj.write_geo_file;
                m = obj.read_msh_file;
                return;
            end

            % get an instance of mesh data base
            m = emdlab_m2d_tmdb;

            % add mesh zones
            for i = 1:numel(obj.faces)
                m.addMeshZone(obj.faces(i).tag, obj.faces(i).getMesh(meshGenerator));
                m.mzs.(obj.faces(i).tag).color = obj.faces(i).color;
            end

        end

        function meshZone = getQMeshByEdges(obj, e1, e2, e3, e4, Nx, Ny)

            % set default values for Nx & Ny
            if nargin < 5, Nx = 3; end
            if nargin < 6, Ny = 3; end

            % set lower limits for Nx & Ny
            Nx = max(Nx+1,3);
            Ny = max(Ny+1,3);

            % edge pointers
            e1ptr = obj.edges(abs(e1)).ptr;
            e2ptr = obj.edges(abs(e2)).ptr;
            e3ptr = obj.edges(abs(e3)).ptr;
            e4ptr = obj.edges(abs(e4)).ptr;

            % Assign nodes to edges
            e1ptr.setNnodes(Nx);
            e2ptr.setNnodes(Ny);
            e3ptr.setNnodes(Nx);
            e4ptr.setNnodes(Ny);

            % Get edge nodes
            pts1 = e1ptr.getMeshNodes;
            pts2 = e2ptr.getMeshNodes;
            pts3 = e3ptr.getMeshNodes;
            pts4 = e4ptr.getMeshNodes;

            % Correct orientation
            if e1<0, pts1 = flipud(pts1); end
            if e2<0, pts2 = flipud(pts2); end
            if e3>0, pts3 = flipud(pts3); end
            if e4>0, pts4 = flipud(pts4); end

            % Parametric grid
            [u,v] = ndgrid(linspace(0,1,Nx),linspace(0,1,Ny));
            u = u(:);
            v = v(:);

            [uIndex,vIndex] = ndgrid(1:Nx,1:Ny);
            uIndex = uIndex(:);
            vIndex = vIndex(:);

            % Corner points
            P00 = pts3(1,:);
            P10 = pts3(end,:);
            P01 = pts1(1,:);
            P11 = pts1(end,:);

            % Coons patch interpolation
            pts = (1-v).*pts3(uIndex,:) + v.*pts1(uIndex,:) + ...
                (1-u).*pts4(vIndex,:) + u.*pts2(vIndex,:) ...
                - ((1-u).*(1-v).*P00 + ...
                u.*(1-v).*P10 + ...
                (1-u).*v.*P01 + ...
                u.*v.*P11);

            % Connectivity list
            index = 0;
            cl = zeros((Nx-1)*(Ny-1),4);

            for j = 1:Ny-1
                for i = 1:Nx-1
                    index = index + 1;

                    n1 = (j-1)*Nx + i;
                    n2 = n1 + 1;
                    n3 = n2 + Nx;
                    n4 = n1 + Nx;

                    cl(index,:) = [n1 n2 n3 n4];
                end
            end

            meshZone = emdlab_m2d_qmz(cl,pts);

        end

        %% visualization methos
        % show the geometry sketch
        function varargout = showSketch(obj, showTags, showWFM)

            if nargin<2
                showTags = true;
                showWFM = false;
            elseif nargin<3
                showWFM = false;
            end

            f = figure('NumberTitle', 'on', 'WindowState', 'maximized', 'name', 'EMDLAB Geometry Visualization', 'color', [0.9,0.9,0.9]);
            hold all;

            % plot points: Np = the number of points
            Np = numel(obj.points);
            p = zeros(Np,2);
            for i = 1:Np
                p(i,1) = obj.points(i).x;
                p(i,2) = obj.points(i).y;
            end
            if ~isempty(p)
                pointTags = cell(1,Np);
                for i = 1:Np
                    pointTags{i} = obj.points(i).tag;
                end
                if showTags
                    text(p(:,1), p(:,2), pointTags, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'BackgroundColor', 'y');
                end
                plot(p(:,1), p(:,2), 's', 'LineWidth', 1.5, 'MarkerEdgeColor','k');
                p_min = min(p,[],1);
                p_max = max(p,[],1);
                plot([min(p_min(1),0),max(p_max(1),0)], [0,0], '--', 'color', [0.2,0.2,0.2]);
                plot([0,0], [min(p_min(2),0),max(p_max(2),0)], '--', 'color', [0.2,0.2,0.2]);
            end

            % plot edges: Ne = the number of edges
            Ne = numel(obj.edges);
            v = cell(Ne,1);
            cl = cell(Ne,1);
            c = zeros(Ne,2);
            for i = 1:Ne
                v{i} = obj.edges(i).ptr.getMeshNodes;
                cl{i} = (1:size(v{i},1)-1)';
                cl{i} = [cl{i},cl{i}+1];
                c(i,:) = obj.edges(i).ptr.getCenter;
            end
            Index = 0;
            for i = 2:Ne
                Index = Index + size(v{i-1},1);
                cl{i} = cl{i} + Index;
            end

            v = cell2mat(v);
            cl = cell2mat(cl);

            if ~isempty(v)
                patch('faces', cl, 'vertices', v, 'edgecolor', 'b', 'linewidth',1.2);
                if showTags

                    edgeTags = cell(1,Ne);
                    for i = 1:Ne
                        edgeTags{i} = obj.edges(i).tag;
                    end
                    text(c(:,1), c(:,2), edgeTags, 'BackgroundColor', 'w', ...
                        'HorizontalAlignment','center','VerticalAlignment','middle');
                end
                if showWFM
                    plot(v(:,1), v(:,2), 'o', 'color', 'k', 'markersize',5, 'markerfacecolor','k');
                end
            end

            set(gca, 'clipping', 'off');

            axis off equal
            zoom on
            grid on
            drawnow;

            if nargout == 1
                varargout{1} = f;
            end

        end

        function varargout = showSketchWithArrows(obj)

f = figure('NumberTitle','on','WindowState','maximized',...
    'name','EMDLAB Geometry Visualization','color',[0.9,0.9,0.9]);
hold all;

%% plot edges (prepare data first)
Ne = numel(obj.edges);
v = cell(Ne,1);
cl = cell(Ne,1);
c = zeros(Ne,2);

for i = 1:Ne
    v{i} = obj.edges(i).ptr.getMeshNodes;
    cl{i} = (1:size(v{i},1)-1)';
    cl{i} = [cl{i},cl{i}+1];
    c(i,:) = obj.edges(i).ptr.getCenter;
end

Index = 0;
for i = 2:Ne
    Index = Index + size(v{i-1},1);
    cl{i} = cl{i} + Index;
end

v = cell2mat(v);
cl = cell2mat(cl);

%% compute directions
if ~isempty(v)

    p1 = v(cl(:,1),:);
    p2 = v(cl(:,2),:);

    midPoints = (p1 + p2)/2;

    dirs = p2 - p1;
    lens = sqrt(sum(dirs.^2,2));
    lens(lens==0) = 1;

    dirs = dirs ./ lens;

    perp = [-dirs(:,2) dirs(:,1)];

    arrowLength = 0.4;
    arrowWidth  = 0.2;

    %% 1️⃣ draw arrows FIRST (bottom layer)
    for i = 1:size(midPoints,1)

        M = midPoints(i,:);
        d = dirs(i,:);
        n = perp(i,:);

        tip = M + arrowLength*d;

        base1 = M + (arrowWidth/2)*n;
        base2 = M - (arrowWidth/2)*n;

        X = [tip(1) base1(1) base2(1)];
        Y = [tip(2) base1(2) base2(2)];

        patch(X,Y,'c','EdgeColor','none');

    end

    %% 2️⃣ draw edges
    patch('faces',cl,'vertices',v,'edgecolor','b','linewidth',1.2);

end


%% plot points
Np = numel(obj.points);
p = zeros(Np,2);

for i = 1:Np
    p(i,1) = obj.points(i).x;
    p(i,2) = obj.points(i).y;
end

if ~isempty(p)

    %% 3️⃣ draw points
    plot(p(:,1),p(:,2),'s','LineWidth',1.5,'MarkerEdgeColor','k');

    p_min = min(p,[],1);
    p_max = max(p,[],1);

    plot([min(p_min(1),0),max(p_max(1),0)],[0,0],'--','color',[0.2,0.2,0.2]);
    plot([0,0],[min(p_min(2),0),max(p_max(2),0)],'--','color',[0.2,0.2,0.2]);

    %% 4️⃣ point tags
    pointTags = cell(1,Np);
    for i = 1:Np
        pointTags{i} = obj.points(i).tag;
    end

    text(p(:,1),p(:,2),pointTags,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top',...
        'BackgroundColor','y');

end


%% 5️⃣ edge tags (top layer)
edgeTags = cell(1,Ne);
for i = 1:Ne
    edgeTags{i} = obj.edges(i).tag;
end

text(c(:,1),c(:,2),edgeTags,...
    'BackgroundColor','w',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle');


set(gca,'clipping','off');

axis off equal
zoom on
grid on
drawnow;

if nargout == 1
    varargout{1} = f;
end

end


        function showFaces(obj)

            m = obj.generateMesh('mm');
            m.showg;

        end
        %% adding primitive loops
        % this function returns loop index and loop handle
        function varargout = addRectangleLoop(obj, x0, y0, w, h)

            p1Index = obj.addPoint(x0,y0);
            p2Index = obj.addPoint(x0+w,y0);
            p3Index = obj.addPoint(x0+w,y0+h);
            p4Index = obj.addPoint(x0,y0+h);

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 0
                obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addCornerRectangleLoop(obj, x0, y0, w, h)

            p1Index = obj.addPoint(x0,y0);
            p2Index = obj.addPoint(x0+w,y0);
            p3Index = obj.addPoint(x0+w,y0+h);
            p4Index = obj.addPoint(x0,y0+h);

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 0
                obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addCenterRectangleLoop(obj, x0, y0, w, h)

            p1Index = obj.addPoint(x0-w/2,y0-h/2);
            p2Index = obj.addPoint(x0+w/2,y0-h/2);
            p3Index = obj.addPoint(x0+w/2,y0+h/2);
            p4Index = obj.addPoint(x0-w/2,y0+h/2);

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 0
                obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = add3PointCornerRectangleLoop(obj, x0, y0, x1, y1, h)

            % check for positive h
            if h < 0
                error('height <h> must be a positive nonzero number.');
            end

            % horizontal vector: width vector
            uw = [x1,y1] - [x0,y0];

            % vertical vector: height vector
            uh = [-uw(2), uw(1)]; uh = h * uh/norm(uh);

            p1Index = obj.addPoint(x0,y0);
            p2Index = obj.addPoint(x1,y1);
            p3Index = obj.addPoint(x1+uh(1),y1+uh(2));
            p4Index = obj.addPoint(x0+uh(1),y0+uh(2));

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 0
                obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = add3PointCenterRectangleLoop(obj, x0, y0, x1, y1, h)

            % check for positive h
            if h < 0
                error('height <h> must be a positive nonzero number.');
            end

            % horizontal vector: width vector
            uw = 2*([x1,y1] - [x0,y0]);

            % vertical vector: height vector
            uh = [-uw(2), uw(1)]; uh = h * uh/norm(uh);

            p1Index = obj.addPoint([x0,y0] - uw/2 - uh/2);
            p2Index = obj.addPoint([x0,y0] + uw/2 - uh/2);
            p3Index = obj.addPoint([x0,y0] + uw/2 + uh/2);
            p4Index = obj.addPoint([x0,y0] - uw/2 + uh/2);

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 0
                obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addParallelogramLoop(obj, x0, y0, w, h)

            p1Index = obj.addPoint(x0,y0);
            p2Index = obj.addPoint(x0+w,y0);
            p3Index = obj.addPoint(x0+w,y0+h);
            p4Index = obj.addPoint(x0,y0+h);

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 0
                obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % this function returns loop index and loop handle
        function varargout = addCircleLoop(obj, x0, y0, r)

            p1Index = obj.addPoint(x0,y0);
            p2Index = obj.addPoint(x0+r,y0);
            p3Index = obj.addPoint(x0-r,y0);

            e1Index = obj.addArc(p1Index, p2Index, p3Index, 1);
            e2Indexe = obj.addArc(p1Index, p3Index, p2Index, 1);

            if nargout == 0
                obj.addLoop(e1Index, e2Indexe);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Indexe);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Indexe);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addAnnularSectorLoop(obj, Ri, Ro, Theta_1, Theta_2, xc, yc)

            % set default center
            if nargin < 6
                xc = 0;
                yc = 0;
            end

            oIndex = obj.addPoint(xc, yc);
            p1Index = obj.addPoint(xc + Ri * cos(Theta_1), yc + Ri * sin(Theta_1));
            p2Index = obj.addPoint(xc + Ro * cos(Theta_1), yc + Ro * sin(Theta_1));
            p3Index = obj.addPoint(xc + Ro * cos(Theta_2), yc + Ro * sin(Theta_2));
            p4Index = obj.addPoint(xc + Ri * cos(Theta_2), yc + Ri * sin(Theta_2));

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addArc(oIndex, p2Index, p3Index, 1);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addArc(oIndex, p4Index, p1Index, 0);

            if nargout == 0
                obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % this function returns loop index and loop handle
        function varargout = addClosedPolylineLoop(obj, x, y)

            Nx = length(x);
            Ny = length(y);

            if Nx~=Ny
                error('Number of x and y coordinates must be the same.');
            end

            p_indices = zeros(1,Nx);
            for i = 1:length(x)
                p_indices(i) = obj.addPoint(x(i),y(i));
            end

            p_indices(end+1) = p_indices(1);
            e_indices = cell(1,Nx);
            for i = 1:Nx
                e_indices{i} = obj.addSegment(p_indices(i),p_indices(i+1));
            end

            if nargout == 0
                obj.addLoop(e_indices{:});
            elseif nargout == 1
                varargout{1} = obj.addLoop(e_indices{:});
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e_indices{:});
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        %% functions to set wireframe mesh
        % setting a fixed mesh length for all geometrical entities
        function setMeshMaxLength(obj, mLength)

            for i = 1:numel(obj.edges)
                obj.edges(i).ptr.setMaxLength(mLength);
            end

            % set mesh size at points for Gmsh
            for i = 1:numel(obj.points)
                obj.points(i).meshSize = mLength;
            end

        end

        % setting a max mesh length for specefied loops
        function setLoopMeshMaxLength(obj, loopIndex, mLength)

            for li = loopIndex

                % get loop handle
                loopHandle = obj.loops(li);

                for i = 1:numel(loopHandle.edges)
                    eptr = loopHandle.edges{i};
                    eptr.setMaxLength(mLength);
                    switch class(eptr)
                        case 'emdlab_g2d_segment'
                            eptr.p0.meshSize = mLength;
                            eptr.p1.meshSize = mLength;
                        case 'emdlab_g2d_arc'
                            eptr.p1.meshSize = mLength;
                            eptr.p2.meshSize = mLength;
                    end
                end

            end

        end

        % setting a max mesh length for specefied edges
        function setEdgeMeshMaxLength(obj, edgeIndex, mLength)

            for ei = edgeIndex

                % get loop handle
                eptr = obj.edges(ei).ptr;
                eptr.setMaxLength(mLength);
                switch class(eptr)
                    case 'emdlab_g2d_segment'
                        eptr.p0.meshSize = mLength;
                        eptr.p1.meshSize = mLength;
                    case 'emdlab_g2d_arc'
                        eptr.p1.meshSize = mLength;
                        eptr.p2.meshSize = mLength;
                end

            end

        end

        % setting a mesh L1 L2 length for specefied edges
        function setEdgeMeshL1L2Length(obj, edgeIndex, L1, L2)

            for ei = edgeIndex

                % get loop handle
                eptr = obj.edges(ei).ptr;
                eptr.setL1L2(L1,L2);
                switch class(eptr)
                    case 'emdlab_g2d_segment'
                        eptr.p0.meshSize = L1;
                        eptr.p1.meshSize = L2;
                    case 'emdlab_g2d_arc'
                        eptr.p1.meshSize = L1;
                        eptr.p2.meshSize = L2;
                end

            end

        end

        % seeting the mesh size using a radial function in cylindrical coordinate system
        function setMeshLengthByRadialFunction(obj, fHandle)

            % check handle function
            if nargin(fHandle) ~= 1
                error('The radial function must get only one input argument.');
            end

            % set mesh size for all edges
            for i = 1:numel(obj.edges)

                if isa(obj.edges(i).ptr, 'emdlab_g2d_segment')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.p0.norm()), fHandle(obj.edges(i).ptr.p1.norm()))
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_arc')
                    obj.edges(i).ptr.setMaxLength(fHandle(obj.edges(i).ptr.p1.norm()));
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_spline')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.pts(1).norm()), fHandle(obj.edges(i).ptr.pts(end).norm()));
                end

            end

            % set mesh size at points for Gmsh
            for i = 1:numel(obj.points)
                obj.points(i).meshSize = fHandle(obj.points(i).norm());
                if isnan(obj.points(i).meshSize), obj.points(i).meshSize = 1; end
            end

        end

        % setting the mesh size using mesh size function f(x) in cartesian coordinate system
        function setMeshLengthByXFunction(obj, fHandle)

            % check handle function
            if nargin(fHandle) ~= 1
                error('The x function must get only one input argument.');
            end

            % set mesh size for all edges
            for i = 1:numel(obj.edges)

                if isa(obj.edges(i).ptr, 'emdlab_g2d_segment')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.p0.x), fHandle(obj.edges(i).ptr.p1.x))
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_arc')
                    obj.edges(i).ptr.setMaxLength(fHandle(obj.edges(i).ptr.p1.x));
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_spline')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.pts(1).x), fHandle(obj.edges(i).ptr.pts(end).x));
                end

            end

            % set mesh size at points for Gmsh
            for i = 1:numel(obj.points)
                obj.points(i).meshSize = fHandle(obj.points(i).x);
                if isnan(obj.points(i).meshSize), obj.points(i).meshSize = 1; end
            end

        end

        % setting the mesh size using mesh size function f(y) in cartesian coordinate system
        function setMeshLengthByYFunction(obj, fHandle)

            % check handle function
            if nargin(fHandle) ~= 1
                error('The y function must get only one input argument.');
            end

            % set mesh size for all edges
            for i = 1:numel(obj.edges)

                if isa(obj.edges(i).ptr, 'emdlab_g2d_segment')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.p0.y), fHandle(obj.edges(i).ptr.p1.y))
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_arc')
                    obj.edges(i).ptr.setMaxLength(fHandle(obj.edges(i).ptr.p1.y));
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_spline')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.pts(1).y), fHandle(obj.edges(i).ptr.pts(end).y));
                end

            end

            % set mesh size at points for Gmsh
            for i = 1:numel(obj.points)
                obj.points(i).meshSize = fHandle(obj.points(i).y);
                if isnan(obj.points(i).meshSize), obj.points(i).meshSize = 1; end
            end

        end

        % setting the mesh size using mesh size function f(x,y) in cartesian coordinate system
        function setMeshLengthByXYFunction(obj, fHandle)

            % check handle function
            if nargin(fHandle) ~= 2
                error('The xy function must get only one input argument.');
            end

            % set mesh size for all edges
            for i = 1:numel(obj.edges)

                % pointer to edge
                eptr = obj.edges(i).ptr;

                if isa(obj.edges(i).ptr, 'emdlab_g2d_segment')
                    obj.edges(i).ptr.setL1L2(fHandle(eptr.p0.x,eptr.p0.y), fHandle(eptr.p1.x,eptr.p1.y))
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_arc')
                    obj.edges(i).ptr.setMaxLength(fHandle(eptr.p1.x,eptr.p1.y));
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_spline')
                    obj.edges(i).ptr.setL1L2(fHandle(eptr.pts(1).x,eptr.pts(1).y), fHandle(eptr.pts(end).x,eptr.pts(end).y));
                end

            end

            % set mesh size at points for Gmsh
            for i = 1:numel(obj.points)
                obj.points(i).meshSize = fHandle(obj.points(i).x,obj.points(i).y);
                if isnan(obj.points(i).meshSize), obj.points(i).meshSize = 1; end
            end

        end

        %% communication with Gmsh software
        function write_geo_file(obj)

            % define a new geo file
            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'w');
            fprintf(fid, 'SetFactory("OpenCASCADE");\n');

            % add points
            for i = 1:numel(obj.points)
                fprintf(fid, 'Point(%d) = {%.16f, %.16f, 0, %f};\n', i, obj.points(i).x, obj.points(i).y, obj.points(i).meshSize);
            end

            % add edges
            for i = 1:numel(obj.edges)

                if isa(obj.edges(i).ptr, 'emdlab_g2d_segment')

                    fprintf(fid, 'Line(%d) = {%d, %d};\n', i, obj.getPointIndexByTag(obj.edges(i).ptr.p0.tag), ...
                        obj.getPointIndexByTag(obj.edges(i).ptr.p1.tag));

                elseif isa(obj.edges(i).ptr, 'emdlab_g2d_arc')

                    fprintf(fid, 'Circle(%d) = {%d, %d, %d};\n', i, obj.getPointIndexByTag(obj.edges(i).ptr.p1.tag), ...
                        obj.getPointIndexByTag(obj.edges(i).ptr.p0.tag), obj.getPointIndexByTag(obj.edges(i).ptr.p2.tag));

                elseif isa(obj.edges(i).ptr, 'emdlab_g2d_spline')

                    pointsList = zeros(1,numel(obj.edges(i).ptr.pts));
                    for j = 1:numel(obj.edges(i).ptr.pts)
                        pointsList(j) = obj.getPointIndexByTag(obj.edges(i).ptr.pts(j).tag);
                    end
                    pointsList = join(string(pointsList), ", ");
                    fprintf(fid, 'Spline(%d) = {%s};\n', i, pointsList);

                end

            end

            % add loops
            for i = 1:numel(obj.loops)

                fprintf(fid, 'Curve Loop(%d) = {%s};\n',i, join(string(obj.loops(i).edgesIndexList), ', '));

            end

            % add faces
            for i = 1:numel(obj.faces)

                tmp_str = num2str(obj.getLoopIndexByTag(obj.faces(i).loops(1).tag));
                for j = 2:length(obj.faces(i).loops)
                    tmp_str = [tmp_str , ', ' , num2str(obj.getLoopIndexByTag(obj.faces(i).loops(j).tag))];
                end

                fprintf(fid, 'Plane Surface(%d) = {%s};\n', i, tmp_str);
                fprintf(fid, 'Physical Surface(%d) = {%d};\n', i, i);

            end

            fclose(fid);

        end

        function m = read_msh_file(obj)

            % run gmsh via matlab
            pyCodePath = "C:\\emdlab-win64\\py-files\\gmsh\\emdlab_gmsh_runGeoSaveMsh2D.py";

            [status ,~] = system(char('"' + obj.pyPath + '"' + " " + '"' + pyCodePath+ '"'));

            % check for system error
            if status ~= 0
                error(['EMDLAB cannot communicate with Gmsh. Please check:\n' ...
                    '1) You have installed Python.\n' ...
                    '2) You have installed Gmsh via: pip install gmsh\n' ...
                    '3) You have set pyPath correctly.\n']);
            end

            % read generated mesh;
            emdlab_gmsh_mshFile;

            % get an instance of mesh data base
            m = emdlab_m2d_tmdb;

            nodes = msh.POS(:,1:2);
            Np = size(nodes,1);

            % add faces
            for i = 1:numel(obj.faces)

                cl = msh.TRIANGLES(msh.TRIANGLES(:,4) == i, 1:3);
                index = unique(cl(:));
                index = sort(index);
                xpoints = nodes(index,:);
                pindex = zeros(Np,1);
                pindex(index) = 1:size(xpoints,1);
                cl = pindex(cl);

                p21 = xpoints(cl(:,2),:) - xpoints(cl(:,1),:);
                p31 = xpoints(cl(:,3),:) - xpoints(cl(:,1),:);
                index = (p21(:,1).*p31(:,2) - p21(:,2).*p31(:,1)) < 0;
                cl(index,:) = cl(index,[1,3,2]);

                m.addMeshZone(obj.faces(i).tag, emdlab_m2d_tmz(cl, xpoints));
                m.mzs.(obj.faces(i).tag).color = obj.faces(i).color;

            end

        end

        function extrudeAndSaveStepSTL(obj, faceName, z1, z2)

            obj.write_geo_file;
            index = obj.getFaceIndexByTag(faceName);

            % define a new geo file
            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'a');
            fprintf(fid, "Extrude {0, 0, %.16f} {Surface{%s};}\n", z2-z1, join(string(1:numel(obj.faces)),','));
            fprintf(fid, "Translate {0, 0, %.16f} {Volume{%d};}\n", z1, index);
            fprintf(fid, "Recursive Delete { Volume{%s};}\n", join(string(setdiff(1:numel(obj.faces), index)),','));
            fprintf(fid, 'Coherence;\n');
            fclose(fid);

            % run gmsh via matlab
            pyCodePath = "C:\\emdlab-win64\\py-files\\emdlab_gmsh_runGeoSaveStep.py";

            [status,~] = system(char('"' + obj.pyPath + '"' + " " + '"' + pyCodePath+ '"'));

            % check for system error
            if status ~= 0
                error(['EMDLAB cannot communicate with Gmsh. Please check:\n' ...
                    '1) You have installed Python.\n' ...
                    '2) You have installed Gmsh via: pip install gmsh\n' ...
                    '3) You have set pyPath correctly.\n']);
            end

            stpPath = "C:\emdlab-win64\geometry\step\emdlab_g3d_stepFile.step";
            copyfile(stpPath, cd + "\" + string(faceName) + ".step")

        end

        %% communication with Maxwell software
        function drawGeometryInMaxwellModel(obj)

            % reference script
            fid1 = fopen('C:\emdlab-win64\geometry\g2d\emdlab_g2d_maxwell.vbs', 'r');

            % modified script
            fid2 = fopen('C:\emdlab-win64\geometry\g2d\emdlab_g2d_maxwellScript.vbs', 'w');

            % read/write loop
            while true

                % check the end of reference file and terminate loop
                if feof(fid1)
                    fclose(fid1);
                    fclose(fid2);
                    break;
                end

                str = fgetl(fid1);

                % detect 'matlab line to import matlab variables
                if strcmpi(str(2:end),'matlab')

                    fprintf(fid2, 'call defineGlobalVariable(oProject, "x_pts", "%s")\n', obj.getPointsXCoordinatesForMaxwell(1:length(obj.points)));
                    fprintf(fid2, 'call makeGBHidden(oProject, "x_pts")\n');
                    fprintf(fid2, 'call defineGlobalVariable(oProject, "y_pts", "%s")\n', obj.getPointsYCoordinatesForMaxwell(1:length(obj.points)));
                    fprintf(fid2, 'call makeGBHidden(oProject, "y_pts")\n');
                    fprintf(fid2, 'call defineGlobalVariable(oProject, "e_angles", "%s")\n', obj.getEdgesAnglesForMaxwell(1:length(obj.edges)));
                    fprintf(fid2, 'call makeGBHidden(oProject, "e_angles")\n');

                    % addfaces
                    for i = 1:numel(obj.faces)

                        faceName = obj.faces(i).tag;
                        lNames = strings(1,numel(obj.faces(i).loops));
                        lIndex = 0;

                        for l = obj.faces(i).loops

                            lIndex = lIndex + 1;
                            lName = faceName + "_loop_" + string(lIndex) + "_";
                            eNames = strings(1,numel(l.edges));

                            % add edges
                            for j = 1:numel(l.edges)

                                eptr = l.edges{j};
                                eNames(j) = lName + eptr.tag;

                                if isa(eptr, 'emdlab_g2d_segment')

                                    fprintf(fid2, 'call drawSegment(oEditor, %d, %d, "%s")\n', obj.getPointIndexByTag(eptr.p0.tag), ...
                                        obj.getPointIndexByTag(eptr.p1.tag), eNames(j));

                                elseif isa(eptr, 'emdlab_g2d_arc')

                                    fprintf(fid2, 'call drawArcCPA(oEditor, %d, %d, %d, "%s")\n', obj.getPointIndexByTag(eptr.p0.tag), ...
                                        obj.getPointIndexByTag(eptr.p1.tag), obj.getEdgeIndexByTag(eptr.tag), eNames(j));

                                elseif isa(eptr, 'emdlab_g2d_spline')

                                    %                             pointsList = zeros(1,numel(obj.edges(i).ptr.pts));
                                    %                             for j = 1:numel(obj.edges(i).ptr.pts)
                                    %                                 pointsList(j) = obj.getPointIndexByTag(obj.edges(i).ptr.pts(j).tag);
                                    %                             end
                                    %                             pointsList = join(string(pointsList), ", ");
                                    %                             fprintf(fid, 'Spline(%d) = {%s};\n', i, pointsList);

                                end

                            end

                            % unite edges
                            fprintf(fid2, 'call uniteEdges(oEditor, "%s")\n', join(eNames,','));

                            % cover loop
                            fprintf(fid2, 'call coverLoop(oEditor, "%s")\n', eNames(1));

                            % save loop name
                            lNames(lIndex) = eNames(1);

                        end

                        % subtract first loop from the rest
                        if numel(lNames) > 1
                            fprintf(fid2, 'call subtract(oEditor, "%s", "%s")\n', join(lNames(2:end),','), lNames(1));
                        end

                        % rename face
                        fprintf(fid2, 'call rename(oEditor, "%s", "%s")\n', lNames(1), faceName);

                        % set face color
                        fprintf(fid2, 'call changeObjectColor(oEditor, "%s", %s)\n', faceName, join(string(floor(obj.faces(i).color*255)),","));

                    end

                else
                    fprintf(fid2, '%s\n', str);
                end

            end

            % run modified script
            system('C:\emdlab-win64\geometry\g2d\emdlab_g2d_maxwellScript.vbs');

        end

        function updateGeometryInMaxwellModel(obj)

            % reference script
            fid1 = fopen('C:\emdlab-win64\geometry\g2d\emdlab_g2d_maxwellUpdate.vbs', 'r');

            % modified script
            fid2 = fopen('C:\emdlab-win64\geometry\g2d\emdlab_g2d_maxwellScript.vbs', 'w');

            % read/write loop
            while true

                % check the end of reference file and terminate loop
                if feof(fid1)
                    fclose(fid1);
                    fclose(fid2);
                    break;
                end

                str = fgetl(fid1);

                % detect 'matlab line to import matlab variables
                if strcmpi(str(2:end),'matlab')

                    fprintf(fid2, 'call updateGlobalVariable(oProject, "x_pts", "%s")\n', obj.getPointsXCoordinatesForMaxwell(1:length(obj.points)));
                    fprintf(fid2, 'call updateGlobalVariable(oProject, "y_pts", "%s")\n', obj.getPointsYCoordinatesForMaxwell(1:length(obj.points)));
                    fprintf(fid2, 'call updateGlobalVariable(oProject, "e_angles", "%s")\n', obj.getEdgesAnglesForMaxwell(1:length(obj.edges)));

                else
                    fprintf(fid2, '%s\n', str);
                end

            end

            % run modified script
            system('C:\emdlab-win64\geometry\g2d\emdlab_g2d_maxwellScript.vbs');

        end

    end

    methods (Static=true)

        %% intersection methods
        % intersection of two infinit lines
        function [xi, yi] = getIntersectionLineLine(x1, y1, ux1, uy1, x2, y2, ux2, uy2)

            % Intersection of two infinite lines defined by:
            % L1: (x1,y1) + t*(ux1,uy1)
            % L2: (x2,y2) + s*(ux2,uy2)

            % Check zero direction vectors
            if (ux1 == 0 && uy1 == 0) || (ux2 == 0 && uy2 == 0)
                error('Zero direction vectors is not acceptable.');
            end

            A = [ux1, -ux2;
                uy1, -uy2];

            b = [x2 - x1;
                y2 - y1];

            % Robust parallel check
            if abs(det(A)) < 1e-12 * (norm([ux1 uy1]) + norm([ux2 uy2]))
                xi = []; yi = [];
                return;
            end

            ts = A \ b;
            t = ts(1);

            % Intersection point on line 1
            xi = x1 + t*ux1;
            yi = y1 + t*uy1;
        end

        % intersection of an infinite line with circle
        function [xi,yi] = getIntersectionLineCircle(x,y,ux,uy,xc,yc,r)
            % Intersection of infinite line:
            %   L(t) = (x, y) + t*(ux, uy)
            % with circle:
            %   (X - xc)^2 + (Y - yc)^2 = r^2

            % Check zero direction vectors
            if abs(ux) < eps && abs(uy) < eps
                error('Zero direction vectors is not acceptable.');
            end

            % Check valid radius
            if r < 0
                error('Circle radius must be a positive double number.');
            end

            % Shift line origin to circle center
            dx = x - xc;
            dy = y - yc;

            % Quadratic coefficients
            A = ux*ux + uy*uy;
            B = 2*(dx*ux + dy*uy);
            C = dx*dx + dy*dy - r*r;

            % Discriminant
            D = B*B - 4*A*C;

            % No real intersection
            if D < 0
                xi = []; yi = [];
                return;
            end

            sqrtD = sqrt(max(D,0));
            t1 = (-B + sqrtD) / (2*A);
            t2 = (-B - sqrtD) / (2*A);

            % Compute intersection points
            x1 = x + t1*ux;  y1 = y + t1*uy;
            x2 = x + t2*ux;  y2 = y + t2*uy;

            % Tangent (one solution)
            if abs(D) < 1e-12 * (A + abs(B) + abs(C))
                xi = x1;
                yi = y1;
                return;
            end

            % Two intersections → return in order of smaller t first
            if abs(t1) <= abs(t2)
                xi = [x1; x2];
                yi = [y1; y2];
            else
                xi = [x2; x1];
                yi = [y2; y1];
            end
        end

        % intersection of an infinite line with an infinite ray
        function [xi, yi] = getIntersectionLineRay(x1, y1, ux1, uy1, x2, y2, ux2, uy2)
            % Intersection of a line and a ray
            % Line: (x1,y1) + t*(ux1, uy1), t in (-inf, inf)
            % Ray:  (x2,y2) + s*(ux2, uy2), s >= 0
            % Returns intersection point (xi, yi) or [] if no intersection on ray

            xi = [];
            yi = [];

            % Check zero direction vectors
            if (ux1 == 0 && uy1 == 0) || (ux2 == 0 && uy2 == 0)
                error('Direction vectors cannot be zero.');
            end

            % Solve for t and s: x1 + t*ux1 = x2 + s*ux2, y1 + t*uy1 = y2 + s*uy2
            A = [ux1, -ux2;
                uy1, -uy2];
            b = [x2 - x1;
                y2 - y1];

            detA = ux1*(-uy2) - uy1*(-ux2);

            % Check if line and ray are parallel
            if abs(detA) < 1e-12
                return; % Parallel, no intersection (or collinear)
            end

            ts = A\b;
            t = ts(1);
            s = ts(2);

            % Intersection must be on ray (s >= 0)
            if s < 0
                return;
            end

            xi = x1 + t*ux1;
            yi = y1 + t*uy1;
        end

        % intersection of an infinite ray with circle
        function [xi,yi] = getIntersectionRayCircle(x,y,ux,uy,xc,yc,r)
            % Intersection of RAY:
            %   L(t) = (x, y) + t*(ux, uy),  t >= 0
            % with circle:
            %   (X - xc)^2 + (Y - yc)^2 = r^2

            % Check zero direction vector
            if abs(ux) < eps && abs(uy) < eps
                error('Zero direction vector is not acceptable.');
            end

            % Check radius
            if r < 0
                error('Circle radius must be a positive number.');
            end

            % Shift coordinate system so circle center is at origin
            dx = x - xc;
            dy = y - yc;

            % Quadratic coefficients for |(dx,dy) + t*(ux,uy)|^2 = r^2
            A = ux*ux + uy*uy;
            B = 2*(dx*ux + dy*uy);
            C = dx*dx + dy*dy - r*r;

            % Discriminant
            D = B*B - 4*A*C;

            % No intersection
            if D < 0
                xi = []; yi = [];
                return;
            end

            % Compute the two potential solutions
            sqrtD = sqrt(max(D,0));
            t1 = (-B + sqrtD) / (2*A);
            t2 = (-B - sqrtD) / (2*A);

            % Accept only t >= 0 (ray condition)
            t_valid = [t1; t2];
            t_valid = t_valid(t_valid >= 0);

            if isempty(t_valid)
                xi = []; yi = [];
                return;
            end

            % Compute intersection points
            xi = x + t_valid * ux;
            yi = y + t_valid * uy;

            % If tangent (one point), keep single output
            if length(t_valid) == 1
                xi = xi(1);
                yi = yi(1);
            end
        end

        % intersection of two circles
        function [xi,yi] = getIntersectionCircleCircle(xc1,yc1,r1,xc2,yc2,r2)

            % Intersection of two circles:
            % (X - xc1)^2 + (Y - yc1)^2 = r1^2
            % (X - xc2)^2 + (Y - yc2)^2 = r2^2

            % Validate radii
            if r1 < 0 || r2 < 0
                error('Circle radii must be positive numbers.');
            end

            % Distance between centers
            dx = xc2 - xc1;
            dy = yc2 - yc1;
            d = sqrt(dx*dx + dy*dy);

            % Check special cases

            % Identical circles → infinite intersections (not solvable uniquely)
            if d < eps && abs(r1 - r2) < eps
                xi = []; yi = [];
                return;
            end

            % No intersection: too far apart or one inside another without touching
            if d > r1 + r2 || d < abs(r1 - r2)
                xi = []; yi = [];
                return;
            end

            % Compute 'a' (distance from circle 1 center to line of intersection)
            a = (r1*r1 - r2*r2 + d*d) / (2*d);

            % Height of intersection points above/below the line between centers
            h_sq = r1*r1 - a*a;
            if h_sq < 0
                h_sq = 0; % numerical safety clamp
            end
            h = sqrt(h_sq);

            % Midpoint between the intersection points
            xm = xc1 + a*dx/d;
            ym = yc1 + a*dy/d;

            % Tangent case → one point
            if abs(h) < 1e-14
                xi = xm;
                yi = ym;
                return;
            end

            % Two intersection points
            rx = -dy * (h/d);
            ry =  dx * (h/d);

            xi = [xm + rx; xm - rx];
            yi = [ym + ry; ym - ry];

        end

        % intersection of two rays
        function [xi, yi] = getIntersectionRayRay(x1, y1, ux1, uy1, x2, y2, ux2, uy2)
            % Intersection of two rays
            % Ray 1: (x1,y1) + t*(ux1, uy1), t >= 0
            % Ray 2: (x2,y2) + s*(ux2, uy2), s >= 0
            % Returns intersection point (xi, yi) or [] if no intersection

            xi = [];
            yi = [];

            % Check zero direction vectors
            if (ux1 == 0 && uy1 == 0) || (ux2 == 0 && uy2 == 0)
                error('Ray direction vectors cannot be zero.');
            end

            % Solve for t and s: x1 + t*ux1 = x2 + s*ux2, y1 + t*uy1 = y2 + s*uy2
            A = [ux1, -ux2;
                uy1, -uy2];
            b = [x2 - x1;
                y2 - y1];

            detA = ux1*(-uy2) - uy1*(-ux2);

            % Check if rays are parallel
            if abs(detA) < 1e-12
                % Parallel rays: no intersection or infinite (collinear)
                % Optional: handle collinear separately if needed
                return;
            end

            ts = A\b;
            t = ts(1);
            s = ts(2);

            % Check if intersection is in forward direction of both rays
            if t < 0 || s < 0
                return; % intersection is behind one of the rays
            end

            xi = x1 + t*ux1;
            yi = y1 + t*uy1;
        end

        % intersection of two finite segments
        function [xi,yi] = getIntersectionSegmentSegment(x1,y1,x2,y2,x3,y3,x4,y4)
            % Returns the intersection point(s) of two finite segments:
            % Segment 1: (x1,y1)-(x2,y2)
            % Segment 2: (x3,y3)-(x4,y4)

            % Direction vectors
            ux = x2 - x1;
            uy = y2 - y1;
            vx = x4 - x3;
            vy = y4 - y3;

            % Solve:
            % (x1,y1) + t*(ux,uy) = (x3,y3) + s*(vx,vy)
            A = [ux, -vx;
                uy, -vy];

            b = [x3 - x1;
                y3 - y1];

            detA = ux*(-vy) - uy*(-vx);

            % Parallel or nearly parallel
            if abs(detA) < 1e-12
                % Check if collinear by checking distance from point to line
                if abs((x3 - x1)*uy - (y3 - y1)*ux) > 1e-12
                    xi = []; yi = [];
                    return; % parallel but not collinear
                end

                % Collinear: check 1D overlap on projection
                % Project onto x or y depending on largest component
                if abs(ux) >= abs(uy)
                    % use x projection
                    seg1 = sort([x1 x2]);
                    seg2 = sort([x3 x4]);

                    left  = max(seg1(1), seg2(1));
                    right = min(seg1(2), seg2(2));

                    if left > right
                        xi = []; yi = [];
                        return; % no overlap
                    end

                    % Overlapping interval in x → compute corresponding points
                    if abs(ux) < 1e-12
                        % vertical line but collinear case handled above
                        xi = x1;
                        yi = linspace(min(y1,y2), max(y1,y2), 2).';
                    else
                        t_left  = (left  - x1) / ux;
                        t_right = (right - x1) / ux;
                        xi = [left; right];
                        yi = [y1 + t_left*uy; y1 + t_right*uy];
                    end

                    return;

                else
                    % use y projection
                    seg1 = sort([y1 y2]);
                    seg2 = sort([y3 y4]);

                    low  = max(seg1(1), seg2(1));
                    high = min(seg1(2), seg2(2));

                    if low > high
                        xi = []; yi = [];
                        return; % no overlap
                    end

                    if abs(uy) < 1e-12
                        yi = y1;
                        xi = linspace(min(x1,x2), max(x1,x2), 2).';
                    else
                        t_low  = (low  - y1) / uy;
                        t_high = (high - y1) / uy;
                        yi = [low; high];
                        xi = [x1 + t_low*ux; x1 + t_high*ux];
                    end

                    return;
                end
            end

            % Non-parallel case → unique intersection if t and s in [0,1]
            ts = A \ b;
            t = ts(1);
            s = ts(2);

            if t < 0 || t > 1 || s < 0 || s > 1
                xi = []; yi = [];
                return; % intersection lies outside segments
            end

            xi = x1 + t*ux;
            yi = y1 + t*uy;
        end

        % intersection of an infinite line with a finite segment
        function [xi,yi] = getIntersectionLineSegment(x,y,ux,uy,x1,y1,x2,y2)
            % Intersection of infinite line:
            %   L(t) = (x, y) + t*(ux, uy)
            % with finite segment:
            %   S(s) = (x1, y1) + s*(dx, dy),  0 <= s <= 1

            % Direction of segment
            dx = x2 - x1;
            dy = y2 - y1;

            % Check zero direction vector for line or segment
            if abs(ux) < eps && abs(uy) < eps
                error('Line direction vector cannot be zero.');
            end

            if abs(dx) < eps && abs(dy) < eps
                error('Segment endpoints are identical; no segment to intersect.');
            end

            % Solve:
            % (x, y) + t*(ux,uy) = (x1, y1) + s*(dx,dy)
            A = [ux, -dx;
                uy, -dy];

            b = [x1 - x;
                y1 - y];

            detA = ux*(-dy) - uy*(-dx);

            % Parallel or nearly parallel
            if abs(detA) < 1e-14
                % Check collinearity
                if abs((x1 - x)*uy - (y1 - y)*ux) > 1e-12
                    xi = []; yi = [];
                    return; % parallel but not collinear
                end

                % Collinear case:
                % Infinite intersections possible → but only return those on segment
                % Parametrize segment onto line: solve s from segment → line param
                % Vector projection of (x1-x,y1-y) onto (ux,uy)
                t1 = ((x1 - x)*ux + (y1 - y)*uy) / (ux*ux + uy*uy);
                t2 = ((x2 - x)*ux + (y2 - y)*uy) / (ux*ux + uy*uy);
                tmin = min(t1, t2);
                tmax = max(t1, t2);

                % Infinite but finite interval → return two endpoints
                xi = [x + tmin*ux; x + tmax*ux];
                yi = [y + tmin*uy; y + tmax*uy];
                return;
            end

            % Non-parallel → unique solution
            ts = A \ b;
            t = ts(1);
            s = ts(2);

            % Check if intersection is on the segment (0 <= s <= 1)
            if s < 0 || s > 1
                xi = []; yi = [];
                return;
            end

            % Compute intersection point
            xi = x + t*ux;
            yi = y + t*uy;
        end

        % intersection of a circle with a finite segment
        function [xi,yi] = getIntersectionCircleSegment(xc, yc, r, x1, y1, x2, y2)
            % Returns intersection points between a circle and a line segment.
            % Circle: center (xc,yc), radius r
            % Segment: endpoints (x1,y1) -> (x2,y2)

            xi = [];
            yi = [];

            % Shift coordinates so circle center becomes origin
            x1s = x1 - xc;
            y1s = y1 - yc;
            x2s = x2 - xc;
            y2s = y2 - yc;

            dx = x2s - x1s;
            dy = y2s - y1s;

            % Quadratic coefficients for intersection with infinite line
            A = dx*dx + dy*dy;
            B = 2*(x1s*dx + y1s*dy);
            C = x1s^2 + y1s^2 - r^2;

            % Discriminant
            D = B*B - 4*A*C;
            if D < 0
                return; % No intersection
            end

            % Compute solutions for t
            sqrtD = sqrt(D);
            t1 = (-B + sqrtD) / (2*A);
            t2 = (-B - sqrtD) / (2*A);

            % Check if each t is within the segment 0 <= t <= 1
            ts = [t1 t2];
            for t = ts
                if t >= 0 && t <= 1
                    xi(end+1) = x1s + t*dx + xc;
                    yi(end+1) = y1s + t*dy + yc;
                end
            end
        end

        % intersection of an infinite ray with a finite segment
        function [xi,yi] = getIntersectionRaySegment(x, y, ux, uy, x1, y1, x2, y2)
            % Intersection of a ray with a line segment.
            % Ray:  start point (x,y), direction (ux,uy) assumed normalized
            % Segment: endpoints (x1,y1) -> (x2,y2)
            %
            % Output:
            %   (xi, yi) intersection point, or [] if no intersection.

            xi = [];
            yi = [];

            % Segment direction
            sx = x2 - x1;
            sy = y2 - y1;

            % Solve for parameters t (ray) and u (segment)
            % Ray:     R(t) = [x; y] + t * [ux; uy],  t >= 0
            % Segment: S(u) = [x1; y1] + u * [sx; sy], 0 <= u <= 1
            %
            % Solve: [x; y] + t[ux;uy] = [x1;y1] + u[sx;sy]

            A = [ux, -sx;
                uy, -sy];

            b = [x1 - x;
                y1 - y];

            detA = A(1,1)*A(2,2) - A(1,2)*A(2,1);
            if abs(detA) < 1e-12
                % Ray and segment are parallel → no intersection
                return;
            end

            % Solve [t; u]
            t = ( b(1)*A(2,2) - A(1,2)*b(2) ) / detA;
            u = ( A(1,1)*b(2) - b(1)*A(2,1) ) / detA;

            % Check ray condition: t >= 0
            if t < 0
                return;
            end

            % Check segment condition: 0 <= u <= 1
            if u < 0 || u > 1
                return;
            end

            % Intersection point
            xi = x + t*ux;
            yi = y + t*uy;

        end

        % intersection of an infinite line with a finite arc
        function [xi, yi] = getIntersectionLineArc(x0, y0, ux, uy, xc, yc, r, theta1, theta2)
            % Intersection of infinite line with a circular arc
            % Line:  (x0,y0) + t*(ux, uy)
            % Arc:   center (xc,yc), radius r, start angle theta1, end angle theta2 (degrees)
            % Returns points lying on the arc

            xi = [];
            yi = [];

            % --- Step 0: check angles ---
            if theta1 >= theta2
                error('theta1 must be less than theta2 (degrees).');
            end

            % Convert degrees to radians
            theta1 = deg2rad(theta1);
            theta2 = deg2rad(theta2);

            % --- Step 1: compute line-circle intersection ---
            dx0 = x0 - xc;
            dy0 = y0 - yc;

            A = ux^2 + uy^2;
            B = 2*(dx0*ux + dy0*uy);
            C = dx0^2 + dy0^2 - r^2;

            D = B^2 - 4*A*C;
            if D < 0
                return; % No intersection
            end

            sqrtD = sqrt(D);
            t1_val = (-B + sqrtD)/(2*A);
            t2_val = (-B - sqrtD)/(2*A);

            % Compute potential intersection points
            points = [x0 + t1_val*ux, y0 + t1_val*uy;
                x0 + t2_val*ux, y0 + t2_val*uy];

            % --- Step 2: filter points to lie on the arc ---
            for k = 1:2
                px = points(k,1);
                py = points(k,2);

                % Compute angle from center to point
                angle = atan2(py - yc, px - xc);

                % Normalize angle to [0,2*pi)
                angle = mod(angle, 2*pi);

                % Check if angle is between theta1 and theta2
                if angle >= theta1 && angle <= theta2
                    xi(end+1,1) = px;
                    yi(end+1,1) = py;
                end
            end
        end

        % intersection of an infinite ray with a finite arc
        function [xi, yi] = getIntersectionRayArc(x0, y0, ux, uy, xc, yc, r, theta1, theta2)
            % Intersection of a ray with a circular arc
            % Ray: (x0,y0) + t*(ux, uy), t >= 0
            % Arc: center (xc,yc), radius r, start/end angles in degrees
            %
            % Returns intersection points that lie on both ray and arc

            xi = [];
            yi = [];

            % --- Step 0: check angles ---
            if theta1 == theta2
                error('theta1 and theta2 must not be equal.');
            end

            % Convert degrees to radians
            theta1 = mod(deg2rad(theta1), 2*pi);
            theta2 = mod(deg2rad(theta2), 2*pi);

            % --- Step 1: compute line-circle intersection ---
            dx0 = x0 - xc;
            dy0 = y0 - yc;

            A = ux^2 + uy^2;
            B = 2*(dx0*ux + dy0*uy);
            C = dx0^2 + dy0^2 - r^2;

            D = B^2 - 4*A*C;
            if D < 0
                return; % No intersection
            end

            sqrtD = sqrt(D);
            t_vals = [(-B + sqrtD)/(2*A), (-B - sqrtD)/(2*A)];

            % --- Step 2: filter intersection points on ray and arc ---
            for t = t_vals
                if t < 0
                    continue; % point is behind the ray
                end

                px = x0 + t*ux;
                py = y0 + t*uy;

                % Compute angle from center to point
                angle = atan2(py - yc, px - xc);
                angle = mod(angle, 2*pi);

                % Check if angle lies on the arc
                if theta1 < theta2
                    on_arc = (angle >= theta1) && (angle <= theta2);
                else
                    % Arc crosses the 2*pi → 0 boundary
                    on_arc = (angle >= theta1) || (angle <= theta2);
                end

                if on_arc
                    xi(end+1,1) = px;
                    yi(end+1,1) = py;
                end
            end
        end

        % intersection of a cicle with a finite arc
        function [xi, yi] = getIntersectionCircleArc(xc1, yc1, r1, xc2, yc2, r2, theta1, theta2)
            % Intersection of a circle and a circular arc
            % Circle: center (xc1,yc1), radius r1
            % Arc: center (xc2,yc2), radius r2, start/end angles in degrees
            %
            % Returns intersection points that lie on the arc

            xi = [];
            yi = [];

            % --- Step 0: check angles ---
            if theta1 == theta2
                error('theta1 and theta2 must not be equal.');
            end

            % Convert degrees to radians
            theta1 = mod(deg2rad(theta1), 2*pi);
            theta2 = mod(deg2rad(theta2), 2*pi);

            % --- Step 1: compute intersection of two circles ---
            dx = xc2 - xc1;
            dy = yc2 - yc1;
            d = sqrt(dx^2 + dy^2);

            % Check for no intersection
            if d > r1 + r2 || d < abs(r1 - r2) || (d==0 && abs(r1-r2)<1e-12)
                return; % no intersection or identical circles
            end

            % Distance from circle1 center to intersection line
            a = (r1^2 - r2^2 + d^2) / (2*d);

            % Height from intersection line to points
            h_sq = r1^2 - a^2;
            h = sqrt(max(h_sq, 0));

            % Midpoint between intersection points
            xm = xc1 + a*dx/d;
            ym = yc1 + a*dy/d;

            % Two intersection points
            rx = -dy * (h/d);
            ry =  dx * (h/d);

            pts = [xm + rx, ym + ry;
                xm - rx, ym - ry];

            % --- Step 2: filter points on the arc ---
            for k = 1:2
                px = pts(k,1);
                py = pts(k,2);

                % Angle from arc center to point
                angle = atan2(py - yc2, px - xc2);
                angle = mod(angle, 2*pi);

                if theta1 < theta2
                    on_arc = (angle >= theta1) && (angle <= theta2);
                else
                    % Arc crosses 2pi → 0 boundary
                    on_arc = (angle >= theta1) || (angle <= theta2);
                end

                if on_arc
                    xi(end+1,1) = px;
                    yi(end+1,1) = py;
                end
            end
        end

        % intersection of a finite segment with a finite arc
        function [xi, yi] = getIntersectionSegmentArc(x1, y1, x2, y2, xc, yc, r, theta1, theta2)
            % Intersection of a segment and a circular arc
            % Segment: (x1,y1) -> (x2,y2)
            % Arc: center (xc,yc), radius r, start/end angles in degrees
            % Returns intersection points lying on both the segment and the arc

            xi = [];
            yi = [];

            % --- Step 0: check angles ---
            if theta1 == theta2
                error('theta1 and theta2 must not be equal.');
            end

            % Convert degrees to radians
            theta1 = mod(deg2rad(theta1), 2*pi);
            theta2 = mod(deg2rad(theta2), 2*pi);

            % --- Step 1: compute intersection of infinite line with circle ---
            dx = x2 - x1;
            dy = y2 - y1;

            % Quadratic coefficients
            x1s = x1 - xc;
            y1s = y1 - yc;

            A = dx^2 + dy^2;
            B = 2*(x1s*dx + y1s*dy);
            C = x1s^2 + y1s^2 - r^2;

            D = B^2 - 4*A*C;
            if D < 0
                return; % no intersection
            end

            sqrtD = sqrt(D);
            t_vals = [(-B + sqrtD)/(2*A), (-B - sqrtD)/(2*A)];

            % --- Step 2: filter points that lie on segment and on arc ---
            for t = t_vals
                if t < 0 || t > 1
                    continue; % outside segment
                end

                px = x1 + t*dx;
                py = y1 + t*dy;

                % Compute angle from arc center to point
                angle = atan2(py - yc, px - xc);
                angle = mod(angle, 2*pi);

                % Check if point is on arc
                if theta1 < theta2
                    on_arc = (angle >= theta1) && (angle <= theta2);
                else
                    % Arc crosses 2pi → 0
                    on_arc = (angle >= theta1) || (angle <= theta2);
                end

                if on_arc
                    xi(end+1,1) = px;
                    yi(end+1,1) = py;
                end
            end
        end

        % intersection of two finite arcs
        function [xi, yi] = getIntersectionArcArc(xc1, yc1, r1, theta11, theta12, xc2, yc2, r2, theta21, theta22)
            % Intersection of two circular arcs
            % Arc1: center (xc1,yc1), radius r1, start/end angles in degrees
            % Arc2: center (xc2,yc2), radius r2, start/end angles in degrees
            % Returns intersection points lying on both arcs

            xi = [];
            yi = [];

            % --- Step 0: check angles ---
            if theta11 == theta12 || theta21 == theta22
                error('Start and end angles must not be equal.');
            end

            % Convert degrees to radians
            theta11 = mod(deg2rad(theta11), 2*pi);
            theta12 = mod(deg2rad(theta12), 2*pi);
            theta21 = mod(deg2rad(theta21), 2*pi);
            theta22 = mod(deg2rad(theta22), 2*pi);

            % --- Step 1: compute circle-circle intersection points ---
            dx = xc2 - xc1;
            dy = yc2 - yc1;
            d = sqrt(dx^2 + dy^2);

            % Check for no intersection
            if d > r1 + r2 || d < abs(r1 - r2) || (d==0 && abs(r1-r2)<1e-12)
                return; % no intersection or identical circles
            end

            % Distance from circle1 center to intersection line
            a = (r1^2 - r2^2 + d^2) / (2*d);

            % Height from line to intersection points
            h_sq = r1^2 - a^2;
            h = sqrt(max(h_sq, 0));

            % Midpoint between intersection points
            xm = xc1 + a*dx/d;
            ym = yc1 + a*dy/d;

            % Two intersection points
            rx = -dy * (h/d);
            ry =  dx * (h/d);

            pts = [xm + rx, ym + ry;
                xm - rx, ym - ry];

            % --- Step 2: filter points on both arcs ---
            for k = 1:2
                px = pts(k,1);
                py = pts(k,2);

                % Angle relative to Arc1 center
                angle1 = atan2(py - yc1, px - xc1);
                angle1 = mod(angle1, 2*pi);

                if theta11 < theta12
                    on_arc1 = (angle1 >= theta11) && (angle1 <= theta12);
                else
                    on_arc1 = (angle1 >= theta11) || (angle1 <= theta12);
                end

                % Angle relative to Arc2 center
                angle2 = atan2(py - yc2, px - xc2);
                angle2 = mod(angle2, 2*pi);

                if theta21 < theta22
                    on_arc2 = (angle2 >= theta21) && (angle2 <= theta22);
                else
                    on_arc2 = (angle2 >= theta21) || (angle2 <= theta22);
                end

                if on_arc1 && on_arc2
                    xi(end+1,1) = px;
                    yi(end+1,1) = py;
                end
            end
        end

        %% point distance from edge objects
        function d = getPointDistanceFromLine(xp, yp, x0, y0, ux, uy)
            u = [ux,uy]; u = u/norm(u);
            p0p = [xp,yp] - [x0,y0];
            d = norm(p0p - dot(p0p,u) * u);
        end

    end

end
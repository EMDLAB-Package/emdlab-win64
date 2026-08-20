% EMDLAB: Electrical Machines Design Laboratory
% emdlab data base class for 2d geometries

classdef emdlab_g2d_db < handle & emdlab_ui_console

    properties

        % points
        points (1,:) emdlab_g2d_point;

        % edges
        edges (1,:) emdlab_g2d_edge;

        % loops
        loops (1,:) emdlab_g2d_loop;

        % faces
        faces (1,:) emdlab_g2d_face;

        % geometrical tolerance
        gtol (1,1) double = 1e-5;

    end

    properties (SetAccess = protected)

        % index holders for ids: identification numbers
        PIDH (1,1) double = 0;
        EIDH (1,1) double = 0;
        LIDH (1,1) double = 0;
        FIDH (1,1) double = 0;

        % id mappers
        pid2pi (1,:) double;
        eid2ei (1,:) double;
        lid2li (1,:) double;
        fid2fi (1,:) double;

        % coordinates for fast check and return
        pts (:,2) double;
        segs (:,2) double;
        arcs (:,2) double;

        % python path
        pyPath = "";

    end

    properties (Dependent = true)
        Npoints (1,1) double;
        Nedges (1,1) double;
        Nloops (1,1) double;
        Nfaces (1,1) double;
    end

    methods
        %% constructor and destructor
        function obj = emdlab_g2d_db()

            obj.printFlag = false;

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

        function y = get.Npoints(obj)
            y = numel(obj.points);
        end

        function y = get.Nedges(obj)
            y = numel(obj.edges);
        end

        function y = get.Nloops(obj)
            y = numel(obj.loops);
        end

        function y = get.Nfaces(obj)
            y = numel(obj.faces);
        end

        function updateIDMappers(obj)

            obj.pid2pi = zeros(1,obj.PIDH);
            obj.eid2ei = zeros(1,obj.EIDH);
            obj.lid2li = zeros(1,obj.LIDH);
            obj.fid2fi = zeros(1,obj.FIDH);
            for i = 1:obj.Npoints
                obj.pid2pi(obj.points(i).id) = i;
            end
            for i = 1:obj.Nedges
                obj.eid2ei(obj.edges(i).id) = i;
            end
            for i = 1:obj.Nloops
                obj.lid2li(obj.loops(i).id) = i;
            end
            for i = 1:obj.Nfaces
                obj.fid2fi(obj.faces(i).id) = i;
            end

        end

        %% point methods
        function varargout = addPoint(obj, varargin)
            % adding a new point to data base
            % this function returns point id and point handle

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
            flg = find(vecnorm(obj.pts - [x,y],2,2) < obj.gtol);

            if flg

                if nargout == 1
                    varargout{1} = obj.points(flg).id;
                elseif nargout == 2
                    varargout{1} = obj.points(flg).id;
                    varargout{2} = obj.points(flg);
                elseif nargout > 2
                    error('The number of output arguments is too high.');
                end
                return;

            end

            % get an instance of point class
            obj.points(end+1) = emdlab_g2d_point(x,y);

            % generate point tag
            obj.PIDH = obj.PIDH + 1;
            obj.points(end).id = obj.PIDH;
            obj.pid2pi(obj.PIDH) = obj.Npoints;

            % add point to database
            obj.pts(end+1,:) = [x,y];

            if nargout == 1
                varargout{1} = obj.PIDH;
            elseif nargout == 2
                varargout{1} = obj.PIDH;
                varargout{2} = obj.points(end);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function removePoints(obj, varargin)

            pIDs = cell2mat(varargin);
            Np = length(pIDs);
            pIndices = zeros(1,Np);

            for i = 1:Np
                pIndices(i) = obj.pid2pi(pIDs(i));
            end

            % remove all edges that are constructed using these points
            elist = cell(1,2*Np);
            for i = 1:Np
                elist{2*i-1} = obj.points(pIndices(i)).ids;
                elist{2*i} = obj.points(pIndices(i)).uids;
            end
            elist = cell2mat(elist);
            obj.removeEdges(elist);

            obj.points(pIndices) = [];
            obj.pts(pIndices,:) = [];
            obj.updateIDMappers;

        end

        function pointIndex = getPointIndexByID(obj, pointID)

            % check for existance of already defined point in data base
            if obj.pid2pi(pointID) == 0
                error('Point with specified ID was not found.');
            else
                pointIndex = obj.pid2pi(pointID);
            end

        end

        function pointHandle = getPointHandleByID(obj, pointID)

            % check for existance of already defined point in data base
            if obj.pid2pi(pointID) == 0
                error('Point with specified ID was not found.');
            else
                pointHandle = obj.points(obj.pid2pi(pointID));
            end

        end

        function pointID = getPointIDByCoordinates(obj, x, y)

            if ~numel(obj.points)
                error('There is no defined point.');
            end

            % this function returns index of closed point to x and y coordinates
            minDistance = inf;
            for i = 1:numel(obj.points)

                distance = norm(obj.points(i).getVector - [x,y]);
                if distance < minDistance
                    pointID = obj.points(i).id;
                    minDistance = distance;
                end

            end

        end

        function [x,y] = getPointCoordinates(obj, pointID)
            pIndex = obj.pid2pi(pointID);
            x = obj.points(pIndex).x;
            y = obj.points(pIndex).y;
        end

        function setPointXCoordinate(obj, pointID, newX)
            obj.points(obj.pid2pi(pointID)).x = newX;
        end

        function setPointYCoordinate(obj, pointID, newY)
            obj.points(obj.pid2pi(pointID)).y = newY;
        end

        function setPointCoordinates(obj, pointID, newX, newY)
            pIndex = obj.pid2pi(pointID);
            obj.points(pIndex).x = newX;
            obj.points(pIndex).y = newY;
        end

        function alignPointsAlongYAxis(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).x = obj.points(obj.pid2pi(varargin{1})).x;
            end

        end

        function alignPointsAlongXAxis(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).y = obj.points(obj.pid2pi(varargin{1})).y;
            end

        end

        function alignPointsAlongRAxis(obj, varargin)

            u_ref = obj.points(obj.pid2pi(varargin{1})).getUnitVector;
            for i = 2:numel(varargin)
                r_i = obj.points(obj.pid2pi(varargin{i})).getDistanceFromOrigin;
                obj.points(obj.pid2pi(varargin{i})).setCoordinates(r_i*u_ref(1), r_i*u_ref(2));
            end

        end

        function alignPointsAlongTAxis(obj, varargin)

            r_ref = obj.points(obj.pid2pi(varargin{1})).getDistanceFromOrigin;
            for i = 2:numel(varargin)
                u_i = obj.points(obj.pid2pi(varargin{i})).getUnitVector;
                obj.points(obj.pid2pi(varargin{i})).setCoordinates(r_ref*u_i(1), r_ref*u_i(2));
            end

        end

        % align points functions
        function alignPointsXX(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).x = obj.points(obj.pid2pi(varargin{1})).x;
            end

        end

        function alignPointsXY(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).y = obj.points(obj.pid2pi(varargin{1})).y;
            end

        end

        function alignPointsXR(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).y = obj.points(obj.pid2pi(varargin{1})).y;
            end

        end

        function alignPointsXT(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).y = obj.points(obj.pid2pi(varargin{1})).y;
            end

        end

        function alignPointsYX(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).y = obj.points(obj.pid2pi(varargin{1})).y;
            end

        end

        function alignPointsYY(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).y = obj.points(obj.pid2pi(varargin{1})).y;
            end

        end

        function alignPointsYR(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).y = obj.points(obj.pid2pi(varargin{1})).y;
            end

        end

        function alignPointsYT(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).y = obj.points(obj.pid2pi(varargin{1})).y;
            end

        end

        function alignPointsRX(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(obj.pid2pi(varargin{i})).y = obj.points(obj.pid2pi(varargin{1})).y;
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

        function u = getUnitVectorP0P1(obj, p0ID, p1ID)
            u = obj.points(obj.pid2pi(p1ID)).getVector - obj.points(obj.pid2pi(p0ID)).getVector;
            u = u/norm(u);
        end

        %% edge methods
        % adding a new segment to data base
        % this function returns edge id and edge handle
        function varargout = addSegment(obj, p0ID, p1ID)

            % check for existance of already defined point in the same location
            for i = 1:obj.Nedges

                if obj.edges(i).isSegment && all(ismember([p0ID, p1ID], obj.edges(i).pid))

                    if nargout == 1
                        varargout{1} = obj.edges(i).id;
                    elseif nargout == 2
                        varargout{1} = obj.edges(i).id;
                        varargout{2} = obj.edges(i);
                    elseif nargout > 2
                        error('The number of output arguments is too high.');
                    end
                    return;

                end

            end

            % get point indices
            p0Index = obj.pid2pi(p0ID);
            p1Index = obj.pid2pi(p1ID);

            % get an instance of the edge class
            edgeHandle = emdlab_g2d_edge;

            % set pointer class to segment
            edgeHandle.ptr = emdlab_g2d_segment(obj.points(p0Index),obj.points(p1Index));
            obj.edges(end+1) = edgeHandle;

            % assign edge id
            obj.EIDH = obj.EIDH + 1;
            edgeHandle.id = obj.EIDH;
            edgeHandle.pid = [p0ID, p1ID];
            obj.eid2ei(obj.EIDH) = obj.Nedges;

            % add edge id to end point ids
            obj.points(p0Index).ids(end+1) = obj.EIDH;
            obj.points(p1Index).ids(end+1) = obj.EIDH;

            if nargout == 1
                varargout{1} = obj.EIDH;
            elseif nargout == 2
                varargout{1} = obj.EIDH;
                varargout{2} = edgeHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % add a new segment by direct coordinates passing
        function varargout = addSegmentByCoordinates(obj, x1, y1, x2, y2)

            p1ID = obj.addPoint(x1,y1);
            p2ID = obj.addPoint(x2,y2);

            if nargout == 0
                obj.addSegment(p1ID,p2ID);
            elseif nargout == 1
                varargout{1} = obj.addSegment(p1ID,p2ID);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegment(p1ID,p2ID);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function addSegmentP0P1XY(obj, p0ID, x, y)

            % add point 1
            p1ID = obj.addPoint(x, y);

            % add segment
            obj.addSegment(p0ID, p1ID);

        end

        function addSegmentP0WH(obj, p0ID, w, h)

            % add point 1
            p0ptr = obj.points(obj.pid2pi(p0ID));
            p1ID = obj.addPoint(p0ptr.x + w, p0ptr.y + h);

            % add segment
            obj.addSegment(p0ID, p1ID);

        end

        % adding a new spline to data base
        % this function returns edge index and edge handle
        function varargout = addSpline(obj, pointIDs)

            % number of points
            Np = length(pointIDs);
            spline_pts = repmat(emdlab_g2d_point, 1, Np);

            for i = 1:Np
                spline_pts(i) = obj.points(obj.pid2pi(pointIDs(i)));
            end

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;

            % set pointer class to segment
            edgeHandle.ptr = emdlab_g2d_spline(spline_pts);
            obj.edges(end+1) = edgeHandle;

            % assign edge id
            obj.EIDH = obj.EIDH + 1;
            edgeHandle.id = obj.EIDH;
            edgeHandle.pid = [spline_pts(1).id, spline_pts(end).id];
            obj.eid2ei(obj.EIDH) = obj.Nedges;

            if nargout == 1
                varargout{1} = obj.EIDH;
            elseif nargout == 2
                varargout{1} = obj.EIDH;
                varargout{2} = edgeHandle;
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
            pointIDs = zeros(1,Nx);
            for i = 1:Nx
                pointIDs(i) = obj.addPoint(x(i),y(i));
            end

            if nargout == 0
                obj.addSpline(pointIDs);
            elseif nargout == 1
                varargout{1} = obj.addSpline(pointIDs);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSpline(pointIDs);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % adding a new arc to data base
        function varargout = addArc(obj, p0ID, p1ID, p2ID, direction)

            % get point indices
            p0Index = obj.pid2pi(p0ID);
            p1Index = obj.pid2pi(p1ID);
            p2Index = obj.pid2pi(p2ID);

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;

            % set pointer class to arc
            if nargin<5, direction = true; end
            edgeHandle.ptr = emdlab_g2d_arc(obj.points(p0Index),obj.points(p1Index),obj.points(p2Index), direction);
            obj.edges(end+1) = edgeHandle;

            % assign edge id
            obj.EIDH = obj.EIDH + 1;
            edgeHandle.id = obj.EIDH;
            edgeHandle.pid = [p1ID, p2ID];
            obj.eid2ei(obj.EIDH) = obj.Nedges;

            % add edge id to end point ids
            obj.points(p0Index).uids(end+1) = edgeHandle.id;
            obj.points(p1Index).ids(end+1) = edgeHandle.id;
            obj.points(p2Index).ids(end+1) = edgeHandle.id;

            if nargout == 1
                varargout{1} = obj.EIDH;
            elseif nargout == 2
                varargout{1} = obj.EIDH;
                varargout{2} = edgeHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addArcByCoordinates(obj, x1, y1, x2, y2, x3, y3, direction)

            p1ID = obj.addPoint(x1,y1);
            p2ID = obj.addPoint(x2,y2);
            p3ID = obj.addPoint(x3,y3);
            if nargin<8, direction = true; end

            if nargout == 0
                obj.addArc(p1ID, p2ID, p3ID, direction);
            elseif nargout == 1
                varargout{1} = obj.addArc(p1ID, p2ID, p3ID, direction);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArc(p1ID, p2ID, p3ID, direction);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addArcByCoordinatesCPA(obj, x1, y1, x2, y2, arcAngle)
            % CPA: center -> point -> arc

            p1ID = obj.addPoint(x1,y1);
            p2ID = obj.addPoint(x2,y2);
            [x3,y3] = emdlab_g2d_rotatePointsXY(x2,y2,arcAngle,x1,y1);
            p3ID = obj.addPoint(x3,y3);
            direction = arcAngle > 0;

            if nargout == 0
                obj.addArc(p1ID, p2ID, p3ID, direction);
            elseif nargout == 1
                varargout{1} = obj.addArc(p1ID, p2ID, p3ID, direction);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArc(p1ID, p2ID, p3ID, direction);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addArcByCoordinatesCPA2(obj, x0, y0, xm, ym, arcAngle)
            % CPA: center -> point -> arc2 -> arc is in both directions

            p0ID = obj.addPoint(x0,y0);
            [x1,y1] = emdlab_g2d_rotatePointsXY(xm,ym,-arcAngle/2,x0,y0);
            p1ID = obj.addPoint(x1,y1);
            [x2,y2] = emdlab_g2d_rotatePointsXY(xm,ym,arcAngle/2,x0,y0);
            p2ID = obj.addPoint(x2,y2);
            direction = arcAngle > 0;

            if nargout == 0
                obj.addArc(p0ID, p1ID, p2ID, direction);
            elseif nargout == 1
                varargout{1} = obj.addArc(p0ID, p1ID, p2ID, direction);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArc(p0ID, p1ID, p2ID, direction);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function addArcP0XYP1A(obj, x0, y0, p1ID, arcAngle)

            p1ptr = obj.points(obj.pid2pi(p1ID));
            obj.addArcByCoordinatesCPA(x0,y0,p1ptr.x,p1ptr.y,arcAngle);

        end

        % get edge handle
        function edgeIndex = getEdgeIndexByID(obj, edgeID)

            % check for existance of already defined edge in data base
            if obj.eid2ei(edgeID) == 0
                error('Edge with specified ID was not found.');
            else
                edgeIndex = obj.eid2ei(edgeID);
            end

        end

        function edgeHandle = getEdgeHandleByID(obj, edgeID)

            % check for existance of already defined edge in data base
            if obj.eid2ei(edgeID) == 0
                error('Edge with specified ID was not found.');
            else
                edgeHandle = obj.edges(obj.eid2ei(edgeID));
            end


        end

        function edgeHandle = getEdgeHandleByTag(obj, eTag)

            % check for existance of already defined edge in data base
            for i = 1:numel(obj.edges)

                if strcmp(obj.edges(i).id,eTag)

                    edgeHandle = obj.edges(i).ptr;
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
                    pIndex = obj.getPointIndexByID(eptr.p0.id);
                case 'emdlab_g2d_arc'
                    pIndex = obj.getPointIndexByID(eptr.p1.id);
                case 'emdlab_g2d_spline'
                    pIndex = obj.getPointIndexByID(eptr.pts(1).id);
            end

        end

        function pIndex = getEdgeEndPointIndex(obj, eIndex)

            % pointer to edge
            eptr = obj.edges(eIndex).ptr;
            switch class(eptr)
                case 'emdlab_g2d_segment'
                    pIndex = obj.getPointIndexByID(eptr.p1.id);
                case 'emdlab_g2d_arc'
                    pIndex = obj.getPointIndexByID(eptr.p2.id);
                case 'emdlab_g2d_spline'
                    pIndex = obj.getPointIndexByID(eptr.pts(end).id);
            end

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

        function [eID, p1ID, p2ID] = connectEdgeCentersBySegment(obj, e1ID, e2ID)
            [~, p1ID] = obj.splitEdge(e1ID);
            [~, p2ID] = obj.splitEdge(e2ID);
            eID = obj.addSegment(p1ID, p2ID);
        end

        function varargout = addParallelSegmentsP0P1(obj, point0ID, point1ID, w)

            % get coordinates
            x0 = obj.pts(obj.pid2pi(point0ID),1);
            y0 = obj.pts(obj.pid2pi(point0ID),2);
            x1 = obj.pts(obj.pid2pi(point1ID),1);
            y1 = obj.pts(obj.pid2pi(point1ID),2);

            % direction & normal vectors
            u = [x1-x0, y1-y0];
            u = u/norm(u);
            v = [-u(2),u(1)];

            p1ID = obj.addPoint(x0-v(1)*w,y0-v(2)*w);
            p2ID = obj.addPoint(x1-v(1)*w,y1-v(2)*w);
            p3ID = obj.addPoint(x0+v(1)*w,y0+v(2)*w);
            p4ID = obj.addPoint(x1+v(1)*w,y1+v(2)*w);

            e1ID = obj.addSegment(p1ID, p2ID);
            e2ID = obj.addSegment(p3ID, p4ID);

            if nargout == 1
                varargout{1} = [e1ID, e2ID, p1ID, p2ID, p3ID, p4ID];
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addParallelSegmentsP0XYP1XY(obj, x0, y0, x1, y1, w)

            % add points
            point0ID = obj.addPoint(x0,y0);
            point1ID = obj.addPoint(x1,y1);

            if nargout == 0
                obj.addParallelSegmentsP0P1(point0ID, point1ID, w);
            elseif nargout == 1
                varargout{1} = obj.addParallelSegmentsP0P1(point0ID, point1ID, w);
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addParallelSegmentsP0P1XY(obj, point0ID, x1, y1, w)

            % add points
            point1ID = obj.addPoint(x1,y1);

            if nargout == 0
                obj.addParallelSegmentsP0P1(point0ID, point1ID, w);
            elseif nargout == 1
                varargout{1} = obj.addParallelSegmentsP0P1(point0ID, point1ID, w);
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addParallelSegmentsP0XYP1(obj, x0, y0, point1ID, w)

            % add points
            point0ID = obj.addPoint(x0,y0);

            if nargout == 0
                obj.addParallelSegmentsP0P1(point0ID, point1ID, w);
            elseif nargout == 1
                varargout{1} = obj.addParallelSegmentsP0P1(point0ID, point1ID, w);
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addRectangleP0P1(obj, point0ID, point1ID)

            p1ID = point0ID;
            p2ID = obj.addPoint(obj.points(obj.pid2pi(point1ID)).x,obj.points(obj.pid2pi(point0ID)).y);
            p3ID = point1ID;
            p4ID = obj.addPoint(obj.points(obj.pid2pi(point0ID)).x,obj.points(obj.pid2pi(point1ID)).y);

            e1ID = obj.addSegment(p1ID, p2ID);
            e2ID = obj.addSegment(p2ID, p3ID);
            e3ID = obj.addSegment(p3ID, p4ID);
            e4ID = obj.addSegment(p4ID, p1ID);

            if nargout == 1
                varargout{1} = [e1ID, e2ID, e3ID, e4ID];
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addRectangleP0P1XY(obj, point0ID, x1, y1)

            % add point 1
            point1ID = obj.addPoint(x1,y1);

            if nargout == 0
                obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout == 1
                varargout{1} = obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addRectangleP0XYP1(obj, x0, y0, point1ID)

            % add point 1
            point0ID = obj.addPoint(x0,y0);

            if nargout == 0
                obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout == 1
                varargout{1} = obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addRectangleP0XYP1XY(obj, x0, y0, x1, y1)

            % add points 0 & 1
            point0ID = obj.addPoint(x0,y0);
            point1ID = obj.addPoint(x1,y1);

            if nargout == 0
                obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout == 1
                varargout{1} = obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addRectangleP0WH(obj, point0ID, W, H)

            % add points 1
            point1ID = obj.addPoint(obj.points(obj.pid2pi(point0ID)).x + W, obj.points(obj.pid2pi(point0ID)).x + H);

            if nargout == 0
                obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout == 1
                varargout{1} = obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addRectangleP0XYWH(obj, x0, y0, W, H)

            % add points 0 & 1
            point0ID = obj.addPoint(x0,y0);
            point1ID = obj.addPoint(x0 + W, y0 + H);

            if nargout == 0
                obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout == 1
                varargout{1} = obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addRectangleP0P1H(obj, x0, y0, W, H)

            % add points 0 & 1
            point0ID = obj.addPoint(x0,y0);
            point1ID = obj.addPoint(x0 + W, y0 + H);

            if nargout == 0
                obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout == 1
                varargout{1} = obj.addRectangleP0P1(point0ID, point1ID);
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addCircle(obj, varargin)

            switch numel(varargin)
                case 2
                    p0ptr = obj.points(obj.pid2pi(varargin{1}));
                    x0 = p0ptr.x;
                    y0 = p0ptr.y;
                    r = varargin{2};

                case 3
                    x0 = varargin{1};
                    y0 = varargin{2};
                    r = varargin{3};

                otherwise
                    error('Wrong number of input arguments.');
            end

            p1Index = obj.addPoint(x0,y0);
            p2Index = obj.addPoint(x0+r,y0);
            p3Index = obj.addPoint(x0-r,y0);

            e1Index = obj.addArc(p1Index, p2Index, p3Index, 1);
            e2Index = obj.addArc(p1Index, p3Index, p2Index, 1);

            if nargout == 1
                varargout{1} = [e1Index, e2Index];
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addCenterRectangle(obj, x0, y0, w, h)

            p1Index = obj.addPoint(x0-w/2,y0-h/2);
            p2Index = obj.addPoint(x0+w/2,y0-h/2);
            p3Index = obj.addPoint(x0+w/2,y0+h/2);
            p4Index = obj.addPoint(x0-w/2,y0+h/2);

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 1
                varargout{1} = [e1Index, e2Index, e3Index, e4Index];
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addPolyline(obj, x, y, closeFlag)

            if nargin <4, closeFlag = false; end

            Nx = length(x);
            Ny = length(y);

            if Nx~=Ny
                error('Number of x and y coordinates must be the same.');
            end

            p_indices = zeros(1,Nx);
            for i = 1:length(x)
                p_indices(i) = obj.addPoint(x(i),y(i));
            end

            if closeFlag
                p_indices(end+1) = p_indices(1);
            else
                Nx = Nx - 1;
            end

            e_indices = zeros(1,Nx);
            for i = 1:Nx
                e_indices(i) = obj.addSegment(p_indices(i),p_indices(i+1));
            end

            if nargout == 1
                varargout{1} = e_indices;
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addAnnularSector(obj, Ri, Ro, Theta_1, Theta_2, xc, yc)

            % set default center
            if nargin < 6
                xc = 0;
                yc = 0;
            end

            oIndex = obj.addPoint(xc, yc);
            p1ID = obj.addPoint(xc + Ri * cos(Theta_1), yc + Ri * sin(Theta_1));
            p2ID = obj.addPoint(xc + Ro * cos(Theta_1), yc + Ro * sin(Theta_1));
            p3ID = obj.addPoint(xc + Ro * cos(Theta_2), yc + Ro * sin(Theta_2));
            p4ID = obj.addPoint(xc + Ri * cos(Theta_2), yc + Ri * sin(Theta_2));

            e1ID = obj.addSegment(p1ID, p2ID);
            e2ID = obj.addArc(oIndex, p2ID, p3ID, 1);
            e3ID = obj.addSegment(p3ID, p4ID);
            e4ID = obj.addArc(oIndex, p4ID, p1ID, 0);

            if nargout == 1
                varargout{1} = [e1ID, e2ID, e3ID, e4ID];
            elseif nargout > 1
                error('The number of output arguments is too high.');
            end

        end

        %% edge edits
        function indicesOfNewEdges = splitEdgeByRatios(obj, edgeIndex, splitRatio)

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
                    p1Index = obj.getPointIndexByID(eptr.p1.id);

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

        function varargout = splitEdge(obj, edgeID, varargin)

            eptr = obj.edges(obj.eid2ei(edgeID));

            if eptr.isSegment

                switch nargout
                    case 0
                        obj.splitSegment(edgeID, varargin{:});
                    case 1
                        varargout{1} = obj.splitSegment(edgeID, varargin{:});
                    case 2
                        [varargout{1}, varargout{2}] = obj.splitSegment(edgeID, varargin{:});
                    otherwise
                        error('Wrong number of output arguments.');
                end

            elseif eptr.isArc

                switch nargout
                    case 0
                        obj.splitArc(edgeID, varargin{:});
                    case 1
                        varargout{1} = obj.splitArc(edgeID, varargin{:});
                    case 2
                        [varargout{1}, varargout{2}] = obj.splitArc(edgeID, varargin{:});
                    otherwise
                        error('Wrong number of output arguments.');
                end

            else
                error('Unsupported edge type.');
            end

        end

        function [newEdgeID, newPointID] = splitSegment(obj, edgeID, pointID)

            % edge & segment pointers
            eptr = obj.edges(obj.eid2ei(edgeID));
            sptr = eptr.ptr;

            if nargin == 2

                tmp = sptr.getCenter;
                newPointID = obj.addPoint(tmp(1),tmp(2));
                p2 = sptr.p1;
                sptr.p1 = obj.points(obj.pid2pi(newPointID));
                newEdgeID = obj.addSegment(newPointID, p2.id);
                obj.points(obj.pid2pi(newPointID)).ids(end+1) = eptr.id;
                p2.ids = setdiff(p2.ids, eptr.id);
                eptr.pid(2) = newPointID;

            elseif nargin == 3

                % point pointer
                pptr = obj.points(obj.pid2pi(pointID));

                if eptr.pid(1) == pptr.id, return; end
                if eptr.pid(2) == pptr.id, return; end
                if ~sptr.isPointOnEdge(pptr); return; end

                %                 p2 = sptr.p1;
                %                 sptr.p1 = pptr;
                %                 newEdgeID = obj.addSegment(pointID, p2.id);
                %                 pptr.ids(end+1) = eptr.id;
                %                 pptr.ids(end+1) = obj.edges(obj.eid2ei(newEdgeID)).id;
                %                 pptr.ids = unique(pptr.ids);
                %                 p2.ids = setdiff(p2.ids, eptr.id);
                %                 eptr.pid(2) = pptr.id;

                pid1 = sptr.p0.id;
                pid2 = sptr.p1.id;

                obj.removeEdges(edgeID);
                newEdgeID(1) = obj.addSegment(pid1, pointID);
                newEdgeID(2) = obj.addSegment(pointID, pid2);

            else
                error('Wrong number of input arguments.');
            end

        end

        function [newEdgeID, newPointID] = splitArc(obj, edgeID, pointID)

            % edge & arc pointers
            eptr = obj.edges(obj.eid2ei(edgeID));
            aptr = eptr.ptr;

            if nargin == 2

                tmp = aptr.p0.getVector;
                tmp = emdlab_g2d_rotatePoints(aptr.p1.getVector, aptr.getSignedAngle/2, tmp(1),tmp(2));
                newPointID = obj.addPoint(tmp(1),tmp(2));
                p2 = aptr.p2;
                aptr.p2 = obj.points(newPointID);
                newEdgeID = obj.addArc(aptr.p0.id, newPointID, p2.id, aptr.direction);
                obj.points(obj.pid2pi(newPointID)).ids(end+1) = eptr.id;
                p2.ids = setdiff(p2.ids, eptr.id);
                eptr.pid(2) = newPointID;

            elseif nargin == 3

                % point pointer
                pptr = obj.points(obj.pid2pi(pointID));

                if eptr.pid(1) == pptr.id, return; end
                if eptr.pid(2) == pptr.id, return; end
                if ~aptr.isPointOnEdge(pptr); return; end

                %                 p2 = aptr.p2;
                %                 aptr.p2 = pptr;
                %                 newEdgeID = obj.addArc(aptr.p0.id, pptr.id, p2.id, aptr.direction);
                %                 pptr.ids(end+1) = eptr.id;
                %                 pptr.ids(end+1) = obj.edges(obj.eid2ei(newEdgeID)).id;
                %                 pptr.ids = unique(pptr.ids);
                %                 p2.ids = setdiff(p2.ids, eptr.id);
                %                 eptr.pid(2) = pptr.id;

                pid0 = aptr.p0.id;
                pid1 = aptr.p1.id;
                pid2 = aptr.p2.id;

                obj.removeEdges(edgeID);
                newEdgeID(1) = obj.addArc(pid0, pid1, pointID, aptr.direction);
                newEdgeID(2) = obj.addArc(pid0, pointID, pid2, aptr.direction);

            else
                error('Wrong number of input arguments.');
            end

        end

        function e = addRectangularFinOnArc(obj, eIndex, wf, hf, t)

            aptr = obj.edges(eIndex).ptr;
            [na1Index, p1Index] = obj.splitArc(eIndex);

            u = obj.points(p1Index).getUnitVector;

            [na2Index, p2Index] = obj.splitArc(na1Index);

            r = aptr.getRadius;
            c = aptr.p0.getVector;

            gamma = 2*asin(wf*0.5/r);

            [x,y] = emdlab_g2d_rotatePointsXY(obj.points(p1Index).x,obj.points(p1Index).y,-gamma/2);
            obj.setPointCoordinates(p1Index, x, y);

            [x,y] = emdlab_g2d_rotatePointsXY(obj.points(p1Index).x,obj.points(p1Index).y,gamma);
            obj.setPointCoordinates(p2Index, x, y);

            p3Index = obj.addPoint(obj.points(p1Index).getVector + hf*u);
            p4Index = obj.addPoint(obj.points(p2Index).getVector + hf*u);

            s1 = obj.addSegment(p1Index, p3Index);
            s2 = obj.addSegment(p3Index, p4Index);
            s3 = obj.addSegment(p4Index, p2Index);

            e = [eIndex, na1Index, na2Index, s1, s2, s3];


        end

        function removeEdges(obj, varargin)

            % id and number of edges that will be removed
            eIDs = cell2mat(varargin);
            Ne = length(eIDs);

            eIndices = zeros(1,Ne);

            for i = 1:Ne

                eIndices(i) = obj.eid2ei(eIDs(i));
                p1 = obj.edges(eIndices(i)).ptr.getPtr1;
                p2 = obj.edges(eIndices(i)).ptr.getPtr2;
                p1.ids = setdiff(p1.ids, eIDs(i));
                p1.uids = setdiff(p1.uids, eIDs(i));
                p2.ids = setdiff(p2.ids, eIDs(i));
                p2.uids = setdiff(p2.uids, eIDs(i));

            end

            obj.edges(eIndices) = [];
            obj.updateIDMappers;

        end

        function removeExcessEdges(obj)

            eIDs = zeros(1,obj.Nedges);
            for i = 1:obj.Nedges
                if isempty(obj.edges(i).ids)
                    eIDs(i) = obj.edges(i).id;
                end
            end
            eIDs(eIDs == 0) = [];
            obj.removeEdges(eIDs);

        end

        function iFlag = intersectEdges(obj, edgeID1, edgeID2)

            iFlag = false;
            [xi,yi] = obj.getIntersection(edgeID1, edgeID2);

            if ~isempty(xi)
                ne = obj.Nedges;
                for i = 1:length(xi)

                    pointID = obj.addPoint(xi(i),yi(i));
                    obj.splitEdge(edgeID1, pointID);
                    obj.splitEdge(edgeID2, pointID);

                    if ne ~= obj.Nedges
                        iFlag = true;
                        break;
                    end
                end
            end

        end

        function intersectAllEdges(obj)

            nfe = 0;
            while true

                existFlag = true;

                for i = 1:obj.Nedges
                    for j = i+1:obj.Nedges
                        nfe = nfe + 1;
                        if obj.intersectEdges(obj.edges(i).id,obj.edges(j).id)
                            existFlag = false;
                        end
                    end
                end

                if existFlag
                    break;
                end

            end

        end

        function buildSketch(obj)

            timeHolder = tic;
            obj.intersectAllEdges;
            obj.dispMessageLine('Intersections completed', timeHolder);

            % remove hanging edges
            timeHolder = tic;
            while true

                elist = [];

                for i = 1:obj.Npoints
                    if length(obj.points(i).ids) == 1
                        elist(end+1) = obj.points(i).ids;
                    end
                end

                if isempty(elist)
                    break;
                end

                obj.removeEdges(elist);

            end
            obj.dispMessageLine('Hanging edges removed', timeHolder);

            % remove unused points
            timeHolder = tic;
            plist = [];

            for i = 1:obj.Npoints
                if isempty(obj.points(i).ids) && isempty(obj.points(i).uids)
                    plist(end+1) = obj.points(i).id;
                end
            end

            obj.removePoints(plist);
            obj.dispMessageLine('Unused points removed', timeHolder);

            obj.updateMMS;

        end

        function updateMMS(obj)
            l_tmp = zeros(1,obj.Nedges);
            for i = 1:obj.Nedges
                if ~obj.edges(i).isSegment
                    l_tmp(i) = obj.edges(i).ptr.getLength();
                else
                    l_tmp(i) = inf;
                end
            end

            obj.setMeshMaxLength(min(l_tmp)/10);
        end

        function applyFillet(obj, fr, pIDs)

            for pID = pIDs

                pptr = obj.points(obj.pid2pi(pID));
                if length(pptr.ids) ~= 2
                    error('To apply fillet on a point, two edges must be connected to it.');
                end
                
                e1ptr = obj.edges(obj.eid2ei(pptr.ids(1)));
                e2ptr = obj.edges(obj.eid2ei(pptr.ids(2)));
                pptr.ids = [];

                eo1ptr = e1ptr.ptr;
                eo2ptr = e2ptr.ptr;

                % check cases
                for shiftValue = [fr,-fr,fr,-fr;fr,fr,-fr,-fr]
                    e1ptr_tmp = eo1ptr.getIndent(shiftValue(1));
                    e2ptr_tmp = eo2ptr.getIndent(shiftValue(2));
                    [xi,yi] = obj.getIntersectionEdges(e1ptr_tmp, e2ptr_tmp);
                    if length(xi) == 1
                        break;
                    end
                end

                if isempty(xi)
                    error('Impossible fillet, fillet radius is too high.');
                end

                p0xy = [xi,yi];
                p1xy = eo1ptr.getPointProjection(xi,yi);
                p2xy = eo2ptr.getPointProjection(xi,yi);

                p0ID = obj.addPoint(p0xy);
                [p1ID, p1PTR] = obj.addPoint(p1xy);
                [p2ID, p2PTR] = obj.addPoint(p2xy);

                u01 = p1xy - p0xy;
                u02 = p2xy - p0xy;

                dirflg = (u01(1) * u02(2) - u01(2) * u02(1)) > 0; 

                obj.addArc(p0ID, p1ID, p2ID, dirflg);

                p1PTR.ids(end+1) = e1ptr.id;
                p2PTR.ids(end+1) = e2ptr.id;                

                if pID == e1ptr.pid(1)
                    eo1ptr.setPtr1(p1PTR);
                    e1ptr.pid(1) = p1ID;
                else
                    eo1ptr.setPtr2(p1PTR);
                    e1ptr.pid(2) = p1ID;
                end

                if pID == e2ptr.pid(1)
                    eo2ptr.setPtr1(p2PTR);
                    e2ptr.pid(1) = p2ID;
                else
                    eo2ptr.setPtr2(p2PTR);
                    e2ptr.pid(2) = p2ID;
                end

            end

        end

        %% loop methods
        % adding a new loop to data base
        % this function returns loop id and loop handle
        function varargout = addLoop(obj, varargin)

            % get a loop class instance
            loopHandle = emdlab_g2d_loop;
            obj.loops(end+1) = loopHandle;

            % loop to assign edges and directions
            for i = 1:numel(varargin)

                for j = 1:numel(varargin{i})

                    edgeID = abs(varargin{i}(j));
                    loopHandle.addEdge(obj.edges(obj.eid2ei(edgeID)).ptr, varargin{i}(j)>0, obj.eid2ei(edgeID));

                end

            end

            % assign loop id
            obj.LIDH = obj.LIDH + 1;
            loopHandle.id = obj.LIDH;
            obj.lid2li(obj.LIDH) = obj.Nloops;

            % add loop id to connected edges
            for edgeIndex = loopHandle.edgesIndexList
                obj.edges(edgeIndex).ids(end+1) = obj.LIDH;
            end

            if nargout == 1
                varargout{1} = obj.LIDH;
            elseif nargout == 2
                varargout{1} = obj.LIDH;
                varargout{2} = loopHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function removeLoop(obj, loopID)

            % get loop index
            loopIndex = obj.lid2li(loopID);

            % first remove all connected faces to this loop
            for faceID = obj.loops(loopIndex).ids
                obj.removeFace(faceID);
            end

            % remove id of this edge from its edges
            for e = obj.loops(loopIndex).edges
                e.ids = setdiff(e.ids, obj.loops(loopIndex).id);
            end

            % remove loop
            obj.loops(loopIndex) = [];

        end

        function removeLoopRecursive(obj, loopIndex)

            if ischar(loopIndex)
                loopIndex = obj.getLoopIndexByTag(loopIndex);
            end

            % first remove all connected faces to this loop
            for faceTag = obj.loops(loopIndex).ids
                obj.removeFace(faceTag);
            end

            obj.loops(loopIndex) = [];

        end

        function loopIndex = getLoopIndexByTag(obj, lTag)

            % check for existance of already defined loop in data base
            for i = 1:numel(obj.loops)

                if strcmp(obj.loops(i).id,lTag)

                    loopIndex = i;
                    return;

                end

            end

            error('Loop was not found.');

        end

        function y = getEdgeLeftLoop(obj, edgeID)
            y = obj.getEdgeLoop(edgeID, true);
        end

        function y = getEdgeRightLoop(obj, edgeID)
            y = obj.getEdgeLoop(edgeID, false);
        end

        function [loopIndex, loopPointer] = getEdgeLoop(obj, edgeID, leftFlag, loopDefineFlag)

            % do you need defining a loop by running this function or not?
            if nargin<4, loopDefineFlag = true; end

            % edge index to start walking
            eIndex = obj.eid2ei(edgeID);

            % get end point of the target edge
            p = obj.points(obj.pid2pi(obj.edges(eIndex).pid(2)));
            elist = eIndex;
            edir = 1;
            pids = obj.edges(eIndex).pid(2);

            % loop for walking
            while true

                % number of edges connected to the point
                n = length(p.ids);
                if n == 1
                    loopIndex = [];
                    loopPointer = [];
                    return;
                end

                % find index of edges connected to the point by ids
                eidx = zeros(1,n);
                for i = 1:n
                    eidx(i) = obj.eid2ei(p.ids(i));
                end

                % remove current edge from edge index list
                eidx = setdiff(eidx, elist(end));
                hflags = false(1,n-1);
                for i = 1:n-1
                    hflags(i) = obj.edges(eidx(i)).isHanging;
                end
                eidx = eidx(~hflags);

                n = length(eidx);
                angles = zeros(1,n);
                flags = zeros(1,n);
                %                 curveture = zeros(1,n);

                eptr = obj.edges(elist(end)).ptr;
                tmp = eptr.getAngles;
                if eptr.isID1(p.id)
                    alpha = tmp(1);
                else
                    alpha = tmp(2);
                end

                %                 if obj.edges(elist(end)).isArc
                %                     crv = 1/eptr.getRadius;
                %                 else
                %                     crv = 0;
                %                 end

                for i = 1:n
                    eptr = obj.edges(eidx(i)).ptr;
                    tmp = eptr.getAngles;
                    %                     if obj.edges(eidx(i)).isArc
                    %                         curveture(i) = 1/eptr.getRadius;
                    %                     else
                    %                         curveture(i) = 0;
                    %                     end

                    if eptr.isID1(p.id)
                        if tmp(1) <= alpha
                            angles(i) = alpha - tmp(1);
                        else
                            angles(i) = 2*pi + alpha - tmp(1);
                        end
                        flags(i) = 1;

                    else
                        if tmp(2) <= alpha
                            angles(i) = alpha - tmp(2);
                        else
                            angles(i) = 2*pi + alpha - tmp(2);
                        end
                        flags(i) = -1;

                    end

                end

                if leftFlag == 1
                    %                     for i = 1:n
                    %                         if abs(angles(i))<1e-5
                    %                             if curveture(i) < crv
                    %                                 angles(i) = 2*pi;
                    %                             end
                    %                         end
                    %                     end
                    %                     tmp = [1:n;angles;-curveture]';
                    %                     tmp = sortrows(tmp, [2,3]);
                    %                     idx = tmp(1,1);
                    [~,idx] = min(angles);
                else
                    %                     for i = 1:n
                    %                         if abs(angles(i))<1e-5
                    %                             if curveture(i) < crv
                    %                                 angles(i) = 2*pi;
                    %                             end
                    %                         end
                    %                     end
                    %                     tmp = [1:n;-angles;curveture]';
                    %                     tmp = sortrows(tmp, [2,3]);
                    %                     idx = tmp(1,1);
                    [~,idx] = max(angles);
                end




                if eidx(idx) == eIndex
                    break;
                end

                if ismember(eidx(idx), elist)
                    loopIndex = [];
                    loopPointer = [];
                    return;
                end

                elist(end+1) = eidx(idx);
                if flags(idx) == 1
                    p = obj.edges(eidx(idx)).ptr.getPtr2;
                    edir(end+1) = 1;
                else
                    p = obj.edges(eidx(idx)).ptr.getPtr1;
                    edir(end+1) = -1;
                end

                %                 if edir(end)>0
                %                 if ismember(obj.edges(eidx(idx)).pid(2),pids)
                % %                             loopIndex = [];
                % %                     loopPointer = [];
                %                     error('Self interseting loops are detected, modify your geometry.');
                % %                     return;
                %                         else
                %                             pids = [pids,obj.edges(eidx(idx)).pid(2)];
                %                 end
                %                 else
                %                     if ismember(obj.edges(eidx(idx)).pid(1),pids)
                %                             %                             loopIndex = [];
                % %                     loopPointer = [];
                %                     error('Self interseting loops are detected, modify your geometry.');
                % %                     return;
                %                         else
                %                             pids = [pids,obj.edges(eidx(idx)).pid(1)];
                %                 end
                %                 end

            end

            if obj.edges(eIndex).ptr.isID2(p.id)
                loopIndex = [];
                loopPointer = [];
                return;
            end

            loopIndex = zeros(1,length(elist));
            for i = 1:length(elist)
                loopIndex(i) = obj.edges(elist(i)).id;
            end

            if ~loopDefineFlag
                loopPointer = [];
                return;
            else
                loopIndex = loopIndex.*edir;
            end

            [loopIndex, loopPointer] = obj.addLoop(loopIndex);

            p = loopPointer.getMeshNodesMinimal;

            pNext = circshift(p, -1, 1);
            signedArea = 0.5 * sum(p(:,1) .* pNext(:,2) - pNext(:,1) .* p(:,2));

            % make it counterclockwise
            if signedArea < 0
                loopPointer.edges = fliplr(loopPointer.edges);
                loopPointer.directions = fliplr(~loopPointer.directions);
                loopPointer.edgesIndexList = fliplr(loopPointer.edgesIndexList);
            end

            % make it canonical
            [~,idx] = min(abs(loopPointer.edgesIndexList));

            loopPointer.edgesIndexList = circshift(loopPointer.edgesIndexList,-idx+1);
            loopPointer.edges = circshift(loopPointer.edges,-idx+1);
            loopPointer.directions = circshift(loopPointer.directions,-idx+1);

        end

        function varargout = updateAllHangingEdges(obj)

            for i = 1:obj.Nedges
                obj.edges(i).isHanging = false;
            end

            hedges = 1:obj.Nedges;
            idx = 1;
            while true

                leftLoop = obj.getEdgeLoop(obj.edges(hedges(idx)).id, true, false);
                rightLoop = obj.getEdgeLoop(obj.edges(hedges(idx)).id, false, false);

                if isempty(leftLoop) && isempty(rightLoop)
                    idx = idx + 1;
                else
                    hedges = setdiff(hedges, [leftLoop, rightLoop]);
                    idx = 1;
                end

                if idx > length(hedges)
                    break;
                end

            end

            for i = hedges
                obj.edges(i).isHanging = true;
            end

            if nargout == 1
                varargout{1} = hedges;
            end

        end

        function removeHangingEdges(obj)
            hedges = obj.updateAllHangingEdges;
            obj.removeEdges(hedges);
        end

        function copyRotateEdges(obj, eIDs, Ncopy, rotAngle)

            if nargin<4, rotAngle = 2*pi/Ncopy; end

            for i = 1:numel(eIDs)

                eptr = obj.edges(obj.eid2ei(eIDs(i)));
                if eptr.isSegment
                    for j = 2:Ncopy
                        [p0x,p0y] = emdlab_g2d_rotatePointsXY(eptr.ptr.p0.x,eptr.ptr.p0.y,(j-1)*rotAngle);
                        [p1x,p1y] = emdlab_g2d_rotatePointsXY(eptr.ptr.p1.x,eptr.ptr.p1.y,(j-1)*rotAngle);
                        obj.addSegmentByCoordinates(p0x,p0y,p1x,p1y);
                    end
                else
                    for j = 2:Ncopy
                        [p0x,p0y] = emdlab_g2d_rotatePointsXY(eptr.ptr.p0.x,eptr.ptr.p0.y,(j-1)*rotAngle);
                        [p1x,p1y] = emdlab_g2d_rotatePointsXY(eptr.ptr.p1.x,eptr.ptr.p1.y,(j-1)*rotAngle);
                        [p2x,p2y] = emdlab_g2d_rotatePointsXY(eptr.ptr.p2.x,eptr.ptr.p2.y,(j-1)*rotAngle);
                        obj.addArcByCoordinates(p0x,p0y,p1x,p1y,p2x,p2y,eptr.ptr.direction);
                    end
                end

            end



        end

        function printLoopEdgeIDs(obj, loopIndex)
            for i = obj.loops(loopIndex).edgesIndexList
                fprintf('%d, ',obj.edges(i).id);
            end
            fprintf('\n');
        end

        %% face methods
        % adding a new face to data base
        % this function returns the face id and face handle
        function varargout = addFace(obj, faceName, varargin)

            % get face class instance
            faceHandle = emdlab_g2d_face;

            % assign face id
            obj.FIDH = obj.FIDH + 1;
            faceHandle.id = obj.FIDH;
            faceHandle.name = faceName;
            faceHandle.color = rand(1,3);
            obj.faces(end+1) = faceHandle;

            for i = 1:numel(varargin)

                for j = 1:numel(varargin{i})
                    loopIndex = obj.lid2li(varargin{i}(j));
                    faceHandle.addLoop(obj.loops(loopIndex));
                    obj.loops(loopIndex).ids(end+1) = obj.FIDH;
                end

            end

            if nargout == 1
                varargout{1} = obj.FIDH;
            elseif nargout == 2
                varargout{1} = obj.FIDH;
                varargout{2} = faceHandle;
            elseif nargout > 2
                error('The number of output arguments is too hight');
            end

        end

        function removeFace(obj, faceID_faceName)

            if ischar(faceID_faceName) || isstring(faceID_faceName)
                faceIndex = obj.getFaceIndexByName(char(faceID_faceName));
            else
                faceIndex = obj.fid2fi(faceID_faceName);
            end

            % remove id of this face from its loops
            for l = obj.faces(faceIndex).loops
                l.ids = setdiff(l.ids, obj.faces(faceIndex).id);
            end

            % remove face
            obj.faces(faceIndex) = [];

        end

        function removeFaceRecursive(obj, faceIndex)

            if ischar(faceIndex)
                faceIndex = obj.getFaceIndexByName(faceIndex);
            end
            for i = 1:obj.faces(faceIndex).Nloops
                obj.removeLoop(obj.faces(faceIndex).loops(i).id);
            end
            obj.faces(faceIndex) = [];

        end

        function faceIndex = getFaceIndexByName(obj, faceName)

            % check for existance of already defined face in data base
            for i = 1:numel(obj.faces)

                if strcmp(obj.faces(i).name, faceName)

                    faceIndex = i;
                    return;

                end

            end

            error('Face was not found.');

        end

        function setFaceColor(obj, faceID_faceName, R, G, B)
            if ischar(faceID_faceName) || isstring(faceID_faceName)
                faceIndex = obj.getFaceIndexByName(char(faceID_faceName));
            else
                faceIndex = obj.fid2fi(faceID_faceName);
            end
            obj.faces(faceIndex).color = [R,G,B]/255;
        end

        function buildFaces(obj)

            timeHolder = tic;
            obj.buildSketch;
            obj.dispMessageLine('Sketch building completed.', timeHolder);

            timeHolder = tic;
            if obj.Nedges == 0
                return;
            end

            idx = 1:obj.Nedges;
            flg = false(1,obj.Nedges);
            for i = 1:obj.Nedges
                flg(i) = obj.edges(i).isHanging;
            end
            idx = idx(~flg);

            l_tmp = cell(1,2*length(idx));
            sli = [];

            rmlist = [];
            for i = 1:length(idx)
                if ismember(idx(i), rmlist)
                    continue
                end
                l_tmp{2*i-1} = obj.getEdgeLoop(obj.edges(idx(i)).id, true, false);
                l_tmp{2*i} = obj.getEdgeLoop(obj.edges(idx(i)).id, false, false);
                if isequal(l_tmp{2*i-1}, l_tmp{2*i}) && ~isempty(l_tmp{2*i})
                    sli(end+1) = idx(i);
                    rmlist = [rmlist, obj.eid2ei(l_tmp{2*i-1})];
                end
            end

            nmax = 0;
            for i = 1:numel(l_tmp)
                nmax = max(nmax, length(l_tmp{i}));
            end

            cl = zeros(numel(l_tmp),nmax);

            for i = 1:numel(l_tmp)
                cl(i,1:length(l_tmp{i})) = l_tmp{i};
            end

            cl(cl(:,1) == 0, :) = [];

            cl = emdlab_mex_makeRowsCanonical(cl);
            %             cl = canonicalizeConnectivity(cl);
            for i = 1:size(cl,1)
                idx = find(cl(i,:) == 0, 1) - 1;
                if idx == 0
                    continue;
                end
                if isempty(idx)
                    idx = size(cl,2);
                end
                if cl(i,2) > cl(i,idx)
                    cl(i,2:idx) = fliplr(cl(i,2:idx));
                end
            end

            cl = unique(cl,'rows');

            tmp = unique(cl(:,1));
            tmp( tmp == 0) = [];

            l_tmpi = [];
            l_tmp = {};

            for i = 1:length(tmp)
                tmp_out = obj.getEdgeLoop(tmp(i), true);
                if ~isempty(tmp_out)
                    l_tmpi(end+1) = tmp_out;
                    l_tmp{end+1} = obj.loops(l_tmpi(end)).edgesIndexList .* (2 * obj.loops(l_tmpi(end)).directions - 1);
                end

                tmp_out =  obj.getEdgeLoop(tmp(i), false);

                if ~isempty(tmp_out)
                    l_tmpi(end+1) = tmp_out;
                    l_tmp{end+1} = obj.loops(l_tmpi(end)).edgesIndexList .* (2 * obj.loops(l_tmpi(end)).directions - 1);
                end

            end

            nmax = 0;
            for i = 1:numel(l_tmp)
                nmax = max(nmax, length(l_tmp{i}));
            end

            cl = zeros(numel(l_tmp),nmax);

            for i = 1:numel(l_tmp)
                cl(i,1:length(l_tmp{i})) = l_tmp{i};
            end

            cl = emdlab_mex_makeRowsCanonical(cl);
            %             cl = canonicalizeConnectivity(cl);
            [~,idx,~] = unique(cl,'rows');


            %
            %             % sort loops vs number of their edges
            %             tmp = zeros(1,length(l_tmpi));
            %             for i = 1:length(l_tmpi)
            %                 tmp(i) = length(l_tmp{i});
            %             end
            %             [~,idx] = sort(tmp);
            %
            %             l_tmpi = l_tmpi(idx);
            %             l_tmp = l_tmp(idx);
            %
            %             % find unique loops
            %             tmp = [];
            %             for i = 1:length(l_tmpi)
            %                 for j = i+1:length(l_tmpi)
            %                     if isequal(l_tmp{i},l_tmp{j})
            %                         tmp(end+1) = i;
            %                         break;
            %                     end
            %                 end
            %             end

            %             idx = setdiff(1:length(l_tmpi),tmp);
            l_tmpi = l_tmpi(idx);
            l_tmp = l_tmp(idx);

            % find boundary loops
            elements = zeros(length(idx),nmax);
            for i = 1:length(idx)
                elements(i,1:length(l_tmp{i})) = l_tmp{i};
            end

            beidx = emdlab_mex_findSignedPairs(elements, max(max(abs(elements))));
            beidx = find(bitor(beidx == 1, beidx == 2));

            % inner face loops
            bidx = [];
            for i = 1:length(idx)
                if all(ismember(abs(l_tmp{i}), beidx))
                    bidx(end+1) = i;
                    continue
                end
            end

            % define special loops
            idx = setdiff(1:length(idx), bidx);

            slidx = [];
            for i = 1:length(bidx)
                if any(ismember(sli,abs(l_tmp{bidx(i)})))
                    slidx(end+1) = bidx(i);
                end
            end

            bidx = setdiff(bidx, slidx);

            % detect boundary loops with zero inner loop
            %             for i = 1:length(bidx)
            %
            %             end

            iloops = {};
            for j = 1:length(idx)
                iloops{j} = l_tmpi(idx(j));
            end

            bloops = {};
            for i = 1:length(bidx)
                bloops{i} = l_tmpi(bidx(i));
            end

            sloops = {};
            for i = 1:length(slidx)
                sloops{i} = l_tmpi(slidx(i));
            end

            for i = 1:length(slidx)
                pts1 = obj.loops(l_tmpi(slidx(i))).getMeshNodesMinimal;
                for j = 1:length(idx)
                    pts2 = obj.loops(l_tmpi(idx(j))).getMeshNodesMinimal;
                    if all(inpolygon(pts1(:,1),pts1(:,2), pts2(:,1),pts2(:,2)))
                        iloops{j}(end+1) = l_tmpi(slidx(i));
                        break;
                    end
                end
            end

            for i = 1:length(bidx)
                pts1 = obj.loops(l_tmpi(bidx(i))).getMeshNodesMinimal;
                for j = 1:length(idx)
                    pts2 = obj.loops(l_tmpi(idx(j))).getMeshNodesMinimal;
                    if all(inpolygon(pts1(:,1),pts1(:,2), pts2(:,1),pts2(:,2)))
                        iloops{j}(end+1) = l_tmpi(bidx(i));
                        break;
                    end
                end
            end

            for i = 1:length(bidx)
                pts1 = obj.loops(l_tmpi(bidx(i))).getMeshNodesMinimal;
                for j = 1:length(slidx)
                    pts2 = obj.loops(l_tmpi(slidx(j))).getMeshNodesMinimal;
                    if all(inpolygon(pts1(:,1),pts1(:,2), pts2(:,1),pts2(:,2)))
                        sloops{j}(end+1) = l_tmpi(bidx(i));
                    end
                end
            end

            for i = 1:length(slidx)
                pts1 = obj.loops(l_tmpi(slidx(i))).getMeshNodesMinimal;
                for j = setdiff(1:length(slidx),i)
                    pts2 = obj.loops(l_tmpi(slidx(j))).getMeshNodesMinimal;
                    if all(inpolygon(pts1(:,1),pts1(:,2), pts2(:,1),pts2(:,2)))
                        sloops{j}(end+1) = l_tmpi(slidx(i));
                    end
                end
            end

            %             for i = 1:length(bidx)
            %                 pts1 = obj.loops(l_tmpi(bidx(i))).getMeshNodesMinimal;
            %                 for j = 1:length(slidx)
            %                     pts2 = obj.loops(l_tmpi(slidx(j))).getMeshNodesMinimal;
            %                     if all(inpolygon(pts1(:,1),pts1(:,2), pts2(:,1),pts2(:,2)))
            %                         sloops{j}(end+1) = l_tmpi(bidx(i));
            %                     end
            %                 end
            %             end

            for i = 1:length(slidx)
                pts1 = obj.loops(l_tmpi(slidx(i))).getMeshNodesMinimal;
                for j = setdiff(1:length(bidx),i)
                    pts2 = obj.loops(l_tmpi(bidx(j))).getMeshNodesMinimal;
                    if all(inpolygon(pts1(:,1),pts1(:,2), pts2(:,1),pts2(:,2)))
                        bloops{j}(end+1) = l_tmpi(slidx(i));
                    end
                end
            end



            allLoops = [iloops, sloops,bloops];

            for i = 1:numel(allLoops)
                if length(allLoops{i}) == 1
                    continue
                end
                childs = [];
                for j = 2:length(allLoops{i})
                    for k = setdiff(1:numel(allLoops),i)
                        if allLoops{k}(1) == allLoops{i}(j)
                            idx = k;
                            break;
                        end
                    end
                    childs = [childs, allLoops{idx}(2:end)];
                end
                allLoops{i} = [allLoops{i}(1), setdiff(allLoops{i}(2:end),childs)];
            end


            fidx = 0;
            for i = 1:(numel(iloops)+numel(sloops))
                fidx = fidx + 1;
                obj.addFace(['f', num2str(fidx)], allLoops{i});
            end

            %             % add faces
            %             fidx = 0;
            %             for i = 1:numel(iloops)
            %                 fidx = fidx + 1;
            %                 obj.addFace(['f', num2str(fidx)], iloops{i});
            %             end
            %
            %             for i = 1:numel(sloops)
            %                 fidx = fidx + 1;
            %                 obj.addFace(['f', num2str(fidx)], sloops{i});
            %             end

            obj.dispMessageLine('Face generation completed.', timeHolder);

        end

        function constructFaces(obj)

            obj.intersectAllEdges;
            obj.updateAllHangingEdges;

            l_tmp = {};
            for i = 1:obj.Nedges
                [~,~,l1] = obj.getEdgeLeftLoop(obj.edges(i).id);
                [~,~,l2] = obj.getEdgeRightLoop(obj.edges(i).id);

                if all(~cellfun(@(x) isequal(x, l1), l_tmp))
                    l_tmp{end+1} = l1;
                end

                if all(~cellfun(@(x) isequal(x, l2), l_tmp))
                    l_tmp{end+1} = l2;
                end

            end

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
                m.addMeshZone(obj.faces(i).name, obj.faces(i).getMesh(meshGenerator));
                m.mzs.(obj.faces(i).name).color = obj.faces(i).color;
            end

        end

        function m = generateQMesh(obj, meshGenerator)

            % default mesh generator
            if nargin<2
                meshGenerator = 'gmsh';
            end

            if strcmpi(meshGenerator, 'gmsh')
                obj.write_geo_file;
                m = obj.read_msh_file_qm;
                return;
            end

            % get an instance of mesh data base
            m = emdlab_m2d_tmdb;

            % add mesh zones
            for i = 1:numel(obj.faces)
                m.addMeshZone(obj.faces(i).id, obj.faces(i).getMesh(meshGenerator));
                m.mzs.(obj.faces(i).id).color = obj.faces(i).color;
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
            meshZone.color = 'c';

        end

        %% visualization methos
        % show the geometry sketch
        function showWireFrameMesh(obj)
            obj.showSketch(0,1);
        end

        function varargout = showSketch(obj, showTagsFlag, showWFMFlag)
            % WFM: wireframe mesh

            if nargin<2, showTagsFlag = true; end
            if nargin<3, showWFMFlag = false; end

            if ~showWFMFlag
                obj.updateMMS;
            end

            figHandle = figure('NumberTitle', 'on', 'name', ...
                'EMDLAB Geometry Visualization', 'color', [0.9,0.9,0.9],'Position',[0,0,1000,600], ...
                'Visible','off');
            movegui(figHandle,'center');
            drawnow;
            figHandle.Visible = 'on';

            figure(figHandle.Number);
            ax = gca;
            cla(ax);
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
                    pointTags{i} = "p" + obj.points(i).id;
                end
                if showTagsFlag
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
                if showTagsFlag

                    edgeTags = cell(1,Ne);
                    for i = 1:Ne
                        edgeTags{i} = "e" + obj.edges(i).id;
                    end
                    text(c(:,1), c(:,2), edgeTags, 'BackgroundColor', 'w', ...
                        'HorizontalAlignment','center','VerticalAlignment','middle');
                end
                if showWFMFlag
                    plot(v(:,1), v(:,2), 'o', 'color', 'k', 'markersize',5, 'markerfacecolor','k');
                end
            end

            set(gca, 'clipping', 'off');

            axis off equal;
            zoom on;
            grid on;
            drawnow;

            if nargout == 1
                varargout{1} = figHandle;
            end

        end

        function varargout = showEdgeDirections(obj)

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

                    patch(X,Y,'r','EdgeColor','none');

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
                    pointTags{i} = "p" + obj.points(i).id;
                end

                text(p(:,1),p(:,2),pointTags,...
                    'HorizontalAlignment','left',...
                    'VerticalAlignment','top',...
                    'BackgroundColor','y');

            end


            %% 5️⃣ edge tags (top layer)
            edgeTags = cell(1,Ne);
            for i = 1:Ne
                edgeTags{i} = "e" + obj.edges(i).id;
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

        function showFaces(obj, varargin)

            m = obj.generateMesh('mm');
            obj.showSketch(0,0);
            ax = gca;
            for i = 1:numel(ax.Children)
                set(ax.Children(i), 'HitTest','off','PickableParts','none');
            end
            m.showg(gca);
            ax = gca;
            ax.Children = flipud(ax.Children);

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
        function setEdgesMeshMaxLength(obj, mLength, varargin)

            eIDs = cell2mat(varargin);

            for eID = eIDs

                ei = obj.eid2ei(eID);

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
        function setEdgesMeshL1L2Length(obj, L1, L2, varargin)

            eIDs = cell2mat(varargin);

            for eID = eIDs

                ei = obj.eid2ei(eID);

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

                    fprintf(fid, 'Line(%d) = {%d, %d};\n', i, obj.getPointIndexByID(obj.edges(i).ptr.p0.id), ...
                        obj.getPointIndexByID(obj.edges(i).ptr.p1.id));

                elseif isa(obj.edges(i).ptr, 'emdlab_g2d_arc')

                    fprintf(fid, 'Circle(%d) = {%d, %d, %d};\n', i, obj.getPointIndexByID(obj.edges(i).ptr.p1.id), ...
                        obj.getPointIndexByID(obj.edges(i).ptr.p0.id), obj.getPointIndexByID(obj.edges(i).ptr.p2.id));

                elseif isa(obj.edges(i).ptr, 'emdlab_g2d_spline')

                    pointsList = zeros(1,numel(obj.edges(i).ptr.pts));
                    for j = 1:numel(obj.edges(i).ptr.pts)
                        pointsList(j) = obj.getPointIndexByID(obj.edges(i).ptr.pts(j).id);
                    end
                    pointsList = join(string(pointsList), ", ");
                    fprintf(fid, 'Spline(%d) = {%s};\n', i, pointsList);

                end

            end

            % add loops
            for i = 1:obj.Nloops

                fprintf(fid, 'Curve Loop(%d) = {%s};\n', i, join(string(obj.loops(i).edgesIndexList), ', '));

            end

            % add faces
            for i = 1:obj.Nfaces

                loopIndexList = zeros(1,obj.faces(i).Nloops);
                for j = 1:length(obj.faces(i).loops)
                    loopIndexList(j) = obj.lid2li(obj.faces(i).loops(j).id);
                end

                fprintf(fid, 'Plane Surface(%d) = {%s};\n', i, strjoin(string(loopIndexList),','));
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

                m.addMeshZone(obj.faces(i).name, emdlab_m2d_tmz(cl, xpoints));
                m.mzs.(obj.faces(i).name).color = obj.faces(i).color;

            end

        end

        function m = read_msh_file_qm(obj)

            % run gmsh via matlab
            pyCodePath = "C:\\emdlab-win64\\py-files\\gmsh\\emdlab_gmsh_runGeoSaveQMsh2D.py";

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
            m = emdlab_m2d_qmdb;

            nodes = msh.POS(:,1:2);
            Np = size(nodes,1);

            % add faces
            for i = 1:numel(obj.faces)

                cl = msh.QUADS(msh.QUADS(:,5) == i, 1:4);
                index = unique(cl(:));
                index = sort(index);
                xpoints = nodes(index,:);
                pindex = zeros(Np,1);
                pindex(index) = 1:size(xpoints,1);
                cl = pindex(cl);

                p21 = xpoints(cl(:,2),:) - xpoints(cl(:,1),:);
                p31 = xpoints(cl(:,3),:) - xpoints(cl(:,1),:);
                index = (p21(:,1).*p31(:,2) - p21(:,2).*p31(:,1)) < 0;
                cl(index,:) = cl(index,[1,4,3,2]);

                m.addMeshZone(obj.faces(i).id, emdlab_m2d_qmz(cl, xpoints));
                m.mzs.(obj.faces(i).id).color = obj.faces(i).color;

            end

        end

        function extrudeAndSaveStepSTL(obj, faceName, z1, z2)

            obj.write_geo_file;
            index = obj.getFaceIndexByName(faceName);

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

                        faceName = obj.faces(i).name;
                        lNames = strings(1,numel(obj.faces(i).loops));
                        lIndex = 0;

                        for l = obj.faces(i).loops

                            lIndex = lIndex + 1;
                            lName = faceName + "_loop_" + string(lIndex) + "_";
                            eNames = strings(1,numel(l.edges));

                            % add edges
                            for j = 1:numel(l.edges)

                                eptr = l.edges{j};
                                eNames(j) = lName + eptr.id;

                                if isa(eptr, 'emdlab_g2d_segment')

                                    fprintf(fid2, 'call drawSegment(oEditor, %d, %d, "%s")\n', obj.getPointIndexByID(eptr.p0.id), ...
                                        obj.getPointIndexByID(eptr.p1.id), eNames(j));

                                elseif isa(eptr, 'emdlab_g2d_arc')

                                    fprintf(fid2, 'call drawArcCPA(oEditor, %d, %d, %d, "%s")\n', obj.getPointIndexByID(eptr.p0.id), ...
                                        obj.getPointIndexByID(eptr.p1.id), obj.getEdgeIndexByID(eptr.id), eNames(j));

                                elseif isa(eptr, 'emdlab_g2d_spline')

                                    %                             pointsList = zeros(1,numel(obj.edges(i).ptr.pts));
                                    %                             for j = 1:numel(obj.edges(i).ptr.pts)
                                    %                                 pointsList(j) = obj.getPointIndexByTag(obj.edges(i).ptr.pts(j).id);
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

        %% get intersection of two edge objects
        function [xi, yi] = getIntersection(obj, e1ID, e2ID)

            e1IDX = obj.eid2ei(e1ID);
            e2IDX = obj.eid2ei(e2ID);
            [xi, yi] = getIntersectionEdges(obj, obj.edges(e1IDX).ptr, obj.edges(e2IDX).ptr);

        end

        function [xi, yi] = getIntersectionEdges(obj, e1ptr, e2ptr)

            switch [class(e1ptr), class(e2ptr)]
                case ['emdlab_g2d_segment', 'emdlab_g2d_segment']
                    [xi,yi] = obj.getIntersectionSegmentSegment(e1ptr.p0.x, e1ptr.p0.y, e1ptr.p1.x, e1ptr.p1.y, ...
                        e2ptr.p0.x, e2ptr.p0.y, e2ptr.p1.x, e2ptr.p1.y);

                case ['emdlab_g2d_segment', 'emdlab_g2d_arc']
                    tmp = e2ptr.getTheta1Theta2;
                    [xi,yi] = obj.getIntersectionSegmentArc(e1ptr.p0.x, e1ptr.p0.y, e1ptr.p1.x, e1ptr.p1.y, ...
                        e2ptr.p0.x, e2ptr.p0.y, e2ptr.getRadius, tmp(1), tmp(2));

                case ['emdlab_g2d_arc', 'emdlab_g2d_segment']
                    tmp = e1ptr.getTheta1Theta2;
                    [xi,yi] = obj.getIntersectionSegmentArc(e2ptr.p0.x, e2ptr.p0.y, e2ptr.p1.x, e2ptr.p1.y, ...
                        e1ptr.p0.x, e1ptr.p0.y, e1ptr.getRadius, tmp(1), tmp(2));

                case ['emdlab_g2d_arc', 'emdlab_g2d_arc']
                    tmp1 = e1ptr.getTheta1Theta2;
                    tmp2 = e2ptr.getTheta1Theta2;
                    [xi,yi] = obj.getIntersectionArcArc(e1ptr.p0.x, e1ptr.p0.y, e1ptr.getRadius, tmp1(1), ...
                        tmp1(2), e2ptr.p0.x, e2ptr.p0.y, e2ptr.getRadius, tmp2(1), tmp2(2));

            end

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
            if abs(h) < 1e-6
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
        function [xi,yi] = getIntersectionSegmentSegment(x1, y1, x2, y2, x3, y3, x4, y4)
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
            if abs(detA) < 1e-5
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
                    if abs(ux) < 1e-5
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

                    if abs(uy) < 1e-5
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
            if d > r1 + r2 || d < abs(r1 - r2) || (d==0 && abs(r1-r2)<1e-5)
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
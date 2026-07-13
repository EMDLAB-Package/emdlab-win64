% EMDLAB: Electrical Machines Design Laboratory
% common properties for all mesh 2d mesh zone database classes

classdef emdlab_m2d_xmdb < handle

    properties

        % mesh nodes: [x,y]
        nodes (:,2) double;

        % mesh connectivity list
        cl (:,:) double;

        % mesh elements: [edge1, edge2, edge3, zone index] -> triangular mesh
        % mesh elements: [edge1, edge2, edge3, edge4, zone index] -> quadrilateral mesh
        elements (:,:) double;

        % unique edges (:,8): [node1, node2, zi1, zi2, ]
        edges

        % list of boundary edges
        bedges

        % edge length
        edgeLength (:,1) double;
        el (:,:) double;
        uEdges (:,2) double;
        nEdges (:,2) double;

        % neighborhood elements
        nbs (:,:) double;

        % jacobian inverse transpose
        JIT (:,:) double;

        % global element area
        gea (1,:) double;

        % element zone index
        ezi (:,:) logical;

        % elements material index
        emi (:,:) logical;

        % auxiliary stored matricies
        mtcs (1,1) struct;

        % named selections
        edgeNamedSelections (1,1) struct;

        % flag to print the elapsed times
        printFlag (1,1) logical = true;

        % element type
        etype (1,:) char = '';

    end

    properties (Dependent = true)

        % Number of nodes
        Nn (1,1) double;

        % Number of elements
        Ne (1,1) double;

    end

    methods

        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end

        function y = get.Ne(obj)
            y = size(obj.cl, 1);
        end

        function setPrintFlag(obj, newValue)
            obj.printFlag = newValue;
        end

        function ggmesh(~)
        end

        %% Visualization Functions
        function varargout = showm(obj, varargin)
            % show global mesh

            [f,ax] = emdlab_flib_fax(varargin{:}); 
            if isfield(f,'MenuBar')
                f.MenuBar = "none";
            end
            
            obj.ggmesh;
            mzNames = string(fieldnames(obj.mzs)');

            for mzName = mzNames
                mzptr = obj.mzs.(mzName);
                if isa(mzptr, 'emdlab_m2d_tmz') || isa(mzptr, 'emdlab_m2d_qmz')
                    plt = patch(ax,'Faces', mzptr.cl, ...
                        'Vertices', mzptr.nodes, 'FaceColor', ...
                        'c', 'EdgeColor', [0.2, 0.2, 0.2], ...
                        'FaceAlpha', 0.7, ...
                        'HitTest','on','PickableParts','visible');
                end
                plt.UserData = mzName;
            end

            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch(ax,'Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5,'HitTest','off','PickableParts','none');

            zoom on;box on;
            ax.Color = [0.86,0.86,0.86];
            grid on;
            grid minor;
            axis(ax, 'equal');
            ax.Toolbar.Visible = 'off';
            ax.GridColor      = [0.7 0.7 0.7];
            ax.MinorGridColor = [0.5 0.5 0.5];
            ax.GridAlpha      = 1;
            ax.MinorGridAlpha = 1;

            set(gcf,'WindowButtonMotionFcn',@hoverFcn);

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

            function hoverFcn(src,~)
                h = hittest(src);
                for i = 1:numel(ax.Children)
                    if isequal(h,ax.Children(i))
                        if isa(ax.Children(i), 'matlab.graphics.primitive.Patch')
                            e.Button = 1;
                            emdlab_flib_selectPatchCallbackGM(ax.Children(i),e);
                            return;
                        end
                    end
                end
                for i = 1:numel(ax.Children)
                    if isa(ax.Children(i), 'matlab.graphics.primitive.Patch')
                        if ischar(ax.Children(i).FaceColor)
                            if strcmpi(ax.Children(i).FaceColor, 'c')
                                set(ax.Children(i), 'FaceColor', 'c', 'FaceAlpha', 0.7);
                                drawnow;
                            end
                        else
                            if any(ax.Children(i).FaceColor ~= [0,1,1])
                                set(ax.Children(i), 'FaceColor', 'c', 'FaceAlpha', 0.7);
                                drawnow;
                            end
                        end
                    end
                end
                title(ax,'');
            end

        end

        function varargout = showgg(obj, varargin)
            % show global geometry

            [f,ax] = emdlab_flib_fax(varargin{:}); 
            if isfield(f,'MenuBar')
                f.MenuBar = "none";
            end
            
            obj.ggmesh;
            mzNames = string(fieldnames(obj.mzs)');

            for mzName = mzNames
                mzptr = obj.mzs.(mzName);
                if isa(mzptr, 'emdlab_m2d_tmz') || isa(mzptr, 'emdlab_m2d_qmz')
                plt = patch(ax,'Faces', mzptr.cl, ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    'c', 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'HitTest','on','PickableParts','visible');
                end
                plt.UserData = mzName;
            end

            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch(ax,'Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5,'HitTest','off','PickableParts','none');

            zoom on;box on;
            ax.Color = [0.86,0.86,0.86];
            grid on;
            grid minor;
            axis(ax, 'equal');
            ax.Toolbar.Visible = 'off';
            ax.GridColor      = [0.4 0.4 0.4];
            ax.MinorGridColor = [0.2 0.2 0.2];
            ax.GridAlpha      = 1;
            ax.MinorGridAlpha = 1;

            set(gcf,'WindowButtonMotionFcn',@hoverFcn);

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

            function hoverFcn(src,~)
                h = hittest(src);
                for i = 1:numel(ax.Children)
                    if isequal(h,ax.Children(i))
                        if isa(ax.Children(i), 'matlab.graphics.primitive.Patch')
                            e.Button = 1;
                            emdlab_flib_selectPatchCallbackGM(ax.Children(i),e);
                            updateXLabel;
                            return;
                        end
                    end
                end
                for i = 1:numel(ax.Children)
                    if isa(ax.Children(i), 'matlab.graphics.primitive.Patch')
                        if ischar(ax.Children(i).FaceColor)
                            if strcmpi(ax.Children(i).FaceColor, 'c')
                                set(ax.Children(i), 'FaceColor', 'c', 'FaceAlpha', 0.5);
                                drawnow;
                            end
                        else
                            if any(ax.Children(i).FaceColor ~= [0,1,1])
                                set(ax.Children(i), 'FaceColor', 'c', 'FaceAlpha', 0.5);
                                drawnow;
                            end
                        end
                    end
                end
                title(ax,'');
                updateXLabel;

                function updateXLabel()
                    % Get cursor position in axes units
                    cp = ax.CurrentPoint;
                    x = cp(1,1);
                    y = cp(1,2);

                    % Check if cursor is inside axes limits
                    xl = ax.XLim;
                    yl = ax.YLim;

                    if x < xl(1) || x > xl(2) || y < yl(1) || y > yl(2)
                        return
                    end

                    % Update xlabel
                    xlabel(ax, sprintf('X = %.2f ,  Y = %.2f ,  R = %.2f ,  D = %.2f',...
                        x, y, norm([x,y]), 2*norm([x,y])), 'Interpreter','none');
                end

            end

        end

        function varargout = showg(obj, varargin)
            % show geometry

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;
            mzNames = fieldnames(obj.mzs);
            if isfield(f,'Name')
                f.Name = "Number of mesh zones = " + string(numel(mzNames));
            end

            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                if isa(mzptr, 'emdlab_m2d_tmz') || isa(mzptr, 'emdlab_m2d_qmz')
                plt = patch('Faces', mzptr.cl, 'Vertices', mzptr.nodes, 'FaceColor', ...
                    mzptr.color, 'EdgeColor', 'none', ...
                    'FaceAlpha', 1, 'Parent', ax);
                end
                plt.UserData.color = mzptr.color;
            end

            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch('Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 0.1, 'parent', ax, 'FaceAlpha', 0.95);

            zoom(ax,'on');
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showfb(obj, varargin)
            % show free boundary

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;

            patch('Faces', obj.edges(obj.bedges, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end
        
        function varargout = showwf(obj, varargin)
            % show wire frame mesh

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;

            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch('Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showmzs(obj, varargin)
            % show mesh zones

            [f,ax] = emdlab_flib_fax(varargin{:});
            mzNames = fieldnames(obj.mzs);
            if isfield(f,'Name')
                f.Name = "Number of mesh zones = " + string(numel(mzNames));
            end

            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                patch('Faces', mzptr.cl, ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    mzptr.color, 'linewidth', 0.05 ,'EdgeColor', [0, 0, 0], ...
                    'FaceAlpha', 1, 'Parent', ax);
            end

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end     
        
        function varargout = showmd(obj, varargin)
            % show mesh degree

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;

            mzNames = fieldnames(obj.mzs);
            ecolor = zeros(obj.Ne, 3);

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                ecolor(obj.ezi(:, i), 1) = mzptr.color(1);
                ecolor(obj.ezi(:, i), 2) = mzptr.color(2);
                ecolor(obj.ezi(:, i), 3) = mzptr.color(3);
            end

            % connectivity list index & point index
            switch obj.etype
                case {'T3', 'TL3'}
                    clIndex = 1:3;
                    pIndex = 1:3;
                case 'TL6'
                    clIndex = 1:3;
                    pIndex = 1:6;
                case 'TL10'
                    clIndex = 1:3;
                    pIndex = 1:10;
                case {'Q4', 'QL4'}
                    clIndex = 1:4;
                    pIndex = 1:4;
                otherwise
                    error('Wrong element type.');
            end

            patch('Faces', obj.cl(:, clIndex), 'Vertices', obj.nodes, ...
                'FaceColor', 'flat', 'FaceVertexCData', ecolor, ...
                'FaceAlpha', 0.5, 'EdgeColor', 'k', 'parent', ax);

            patch('Faces', obj.cl(:, pIndex), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'none', 'parent', ax, 'Marker', 'o', 'MarkerFaceColor', 'k');

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showNodes(obj, varargin)
            % show nodes on global mesh

            [f,ax] = obj.showm(varargin{1});
            if isnumeric(varargin{1})
                sIndex = 1;
            else
                sIndex = 2; 
            end
            color = 'r';
            for i = sIndex:numel(varargin)
                patch('Faces', varargin{i}, 'Vertices', obj.nodes, ...
                    'EdgeColor', 'none', 'parent', ax, 'Marker', 'o', 'MarkerFaceColor', color);
                color = 'b';
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showEdges(obj, varargin)
            % show edges on global mesh

            [f,ax] = obj.showm;
            for i = 1:numel(varargin)
                patch('Faces', obj.edges(varargin{i},[1,2]), 'Vertices', obj.nodes, ...
                    'EdgeColor', 'r', 'parent', ax, 'Marker', 'o', 'MarkerFaceColor', 'r', 'LineWidth', 2);
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showElements(obj, varargin)
            % show elements on global mesh

            [f,ax] = obj.showm;

            % connectivity list index & point index
            switch obj.etype
                case {'T3', 'TL3'}
                    clIndex = 1:3;
                case 'TL6'
                    clIndex = 1:3;
                case 'TL10'
                    clIndex = 1:3;
                case {'Q4', 'QL4'}
                    clIndex = 1:4;
                otherwise
                    error('Wrong element type.');
            end

            for i = 1:numel(varargin)
                patch('Faces', obj.cl(varargin{i}, clIndex), 'Vertices', obj.nodes, ...
                    'FaceColor', 'r', 'EdgeColor', 'r', 'parent', ax, 'Marker', 'o', ...
                    'MarkerFaceColor', 'r', 'LineWidth', 2,'HitTest','off', 'PickableParts', 'none');
            end

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = showce(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            obj.ggmesh;

            pts = [obj.getCenterOfElements; obj.getCenterOfEdges];
            cl1_tmp = repmat(1:obj.Ne,3,1);
            cl2_tmp = abs(obj.elements(:,1:3))' + obj.Ne;
            cl3_tmp = [cl1_tmp(:),cl2_tmp(:)];

            patch('Faces', obj.cl(:, [1,2,3]), 'Vertices', obj.nodes, ...
                'FaceColor', 'w', ...
                'FaceAlpha', 1, 'EdgeColor', 'k', 'linewidth', 1, 'parent', ax);

            patch('Faces', cl3_tmp, 'Vertices', pts, ...
                'FaceColor', 'none', 'EdgeColor', 'r', 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = shownes(obj)
            f = GraphicWindow;
            f.Name = '[Named edges of mesh]';
            h = guihandles(f);
            nes = fieldnames(obj.edgeNamedSelections);

            for i = 1:numel(nes)
                tmp = rand(1, 3);
                patch('Faces', obj.edges(obj.edgeNamedSelections.(nes{i}), 1:2), 'Vertices', ...
                    obj.nodes, 'FaceColor', ...
                    tmp, 'EdgeColor', tmp, 'linewidth', 1.5, ...
                    'FaceAlpha', 1, 'parent', h.va);
            end

            legend(h.va, nes);
            set(f, 'HandleVisibility', 'off', 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function showne(obj, neName)
            neName = obj.checkEdgeNamedSelectionExistence(neName);
            ah = setFigure('TMDBC: Named Edge');
            patch('Faces', obj.edges(obj.edgeNamedSelections.(neName), 1:2), 'Vertices', ...
                obj.nodes, 'FaceColor', ...
                'b', 'EdgeColor', 'b', 'linewidth', 1.5, ...
                'FaceAlpha', 1, 'parent', ah);
            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end

        function varargout = showvf(obj, Fx, Fy)
            f = GraphicWindow;
            f.Name = 'Field plot on center of elements';
            h = guihandles(f);
            h.va.NextPlot = 'add';
            c = obj.getCenterOfElements;
            color = zeros(obj.Ne, 3);
            mzNames = fieldnames(obj.mzs);

            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                color(obj.ezi(:, mzptr.zi), :) = repmat(mzptr.color, mzptr.Ne, 1);
            end

            patch(h.va, 'faces', obj.cl(:, 1:3), 'vertices', obj.nodes, ...
                'facecolor', 'flat', 'FaceVertexCData', color, 'EdgeColor', 'w', 'facealpha', 0.5);
            quiver(h.va, c(:, 1), c(:, 2), Fx, Fy, 'color', 'k');
            axis(h.va, 'off', 'equal');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end

        function varargout = showContact(obj, mz1Name, mz2Name)

            obj.ggmesh;
            idx = obj.getEdgeIndicesOnContact(mz1Name, mz2Name);
            [f,ax] = obj.showgg;

            patch('Faces', obj.edges(idx,1:2), 'Vertices', obj.nodes, ...
                'EdgeColor', 'b', 'parent', ax, 'Marker', 'none', 'MarkerFaceColor', 'none', ...
                'LineWidth', 3, 'PickableParts','none');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        %% Copy & Transform Functions
        function mirrorMeshZone(obj, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(mzName).mirror(varargin{:});
        end

        function mmz(varargin)
            mirrorMeshZone(varargin{:});
        end

        function copyMirrorMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getMirror(varargin{:});
        end

        function cmmz(varargin)
            copyMirrorMeshZone(varargin{:});
        end

        function rotateMeshZone(obj, mzName, varargin)
            mzName = char(mzName);
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(mzName).rotate(varargin{:});
        end

        function rotateMeshZones(obj, mzNames, varargin)
            for mzName = mzNames
                obj.rotateMeshZone(mzName, varargin{:});
            end
        end

        function rmz(varargin)
            rotateMeshZone(varargin{:});
        end

        function copyRotateMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getRotate(varargin{:});
        end

        function crmz(varargin)
            copyRotateMeshZone(varargin{:});
        end

        function shiftMeshZone(obj, mzName, varargin)
            mzName = erase(mzName, ' ');
            obj.checkMeshZoneExistence(mzName)
            obj.mzs.(mzName).shift(varargin{:});
        end

        function shiftMeshZones(obj, mzNames, shiftX, shiftY)
            for mzName = mzNames
                obj.shiftMeshZone(mzName, [shiftX, shiftY]);
            end
        end

        function shmz(varargin)
            shiftMeshZone(varargin{:});
        end

        function copyShiftMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getShift(varargin{:});
        end

        function cshmz(varargin)
            copyShiftMeshZone(varargin{:});
        end

        %% Boolean Functions
        function joinMeshZones(obj, nmzName, varargin)

            % find total number of mesh zones need to be joined
            xNmzs = 0;
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    xNmzs = xNmzs + 1;
                elseif isstring(varargin{i})
                    xNmzs = xNmzs + numel(varargin{i});
                else
                    error('Input type must be <char> or <string>.');
                end
            end

            if xNmzs < 2
                error('Minimum number mzs must be 2.');
            end

            mzNames = cell(1,xNmzs);
            index = 0;
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    index = index + 1;
                    mzNames{index} = obj.checkMeshZoneExistence(varargin{i});
                else
                    for j = 1:numel(varargin{i})
                        index = index + 1;
                        mzNames{index} = obj.checkMeshZoneExistence(char(varargin{i}(j)));
                    end
                end
            end

            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            Nn_tmp = zeros(1, xNmzs);
            Ne_tmp = zeros(1, xNmzs);

            for i = 1:xNmzs
                Nn_tmp(i) = obj.mzs.(mzNames{i}).Nn;
                Ne_tmp(i) = obj.mzs.(mzNames{i}).Ne;
            end

            n_nmz = zeros(sum(Nn_tmp), 2);
            e_nmz = zeros(sum(Ne_tmp), 4);
            n_tmp = 0;
            e_tmp = 0;

            for i = 1:xNmzs
                n_nmz(1 + n_tmp:n_tmp + Nn_tmp(i), :) = obj.mzs.(mzNames{i}).nodes;
                e_nmz(1 + e_tmp:e_tmp + Ne_tmp(i), :) = obj.mzs.(mzNames{i}).cl + n_tmp;
                n_tmp = n_tmp + Nn_tmp(i);
                e_tmp = e_tmp + Ne_tmp(i);
            end

            % jointing mzs
            [n_nmz, ~, ic] = uniquetol(n_nmz, obj.gleps, 'ByRows', true);
            e_nmz = ic(e_nmz);
            % adding new mz
            obj.mzs.(nmzName) = emdlab_m2d_qmz(e_nmz, n_nmz);
            obj.mzs.(nmzName).material = obj.mzs.(mzNames{1}).material;
            obj.mzs.(nmzName).color = obj.mzs.(mzNames{1}).color;
            % removing old mzs
            for i = 1:xNmzs
                obj.mzs = rmfield(obj.mzs, mzNames{i});
            end

        end

        function jmzs(varargin)
            joinMeshZones(varargin{:});
        end

        %% Index Finding Functions
        function y = getfbe(obj)
            % get indices of free boundary edges
            obj.ggmesh;
            y = find(obj.bedges);
        end

        function idx = getNodeIndicesOnLineP0P1(obj, x0, y0, x1, y1, tol)
            % getNodeIndicesOnLineP0P1
            % Returns node indices located within tol of the line segment P0→P1.

            % Handle default tolerance
            if nargin < 6
                tol = obj.gleps;
            end

            % Handle default P1 (if user only gives one point)
            if nargin < 5
                x1 = 0;
                y1 = 0;
            end

            % Construct points
            p1 = [x0, y0];
            p2 = [x1, y1];

            % Direction unit vector
            d = p2 - p1;
            L = norm(d);

            if L < eps
                error('P0 and P1 are identical – line definition invalid.');
            end

            u = d / L;       % unit direction

            % Node coordinates relative to P1
            rel = obj.nodes(:,1:2) - p1;

            % Projection parameter alpha
            alpha = rel * u';     % dot product with direction

            % Perpendicular distance from line
            dist = sqrt(sum((rel - alpha*u).^2, 2));

            % Indices of nodes near the line
            idx = find(dist < tol);
        end

        function idx = getNodeIndicesOnLineP0U(obj, x0, y0, ux, uy, tol)
            % getNodeIndicesOnLineP0U
            % Returns node indices near the line passing through P0 = (x0,y0)
            % in direction U = (ux,uy)

            % Default tolerance
            if nargin < 6
                tol = obj.gleps;
            end

            % Check direction vector
            if ux == 0 && uy == 0
                error('Direction vector U = (ux, uy) must be nonzero.');
            end

            % Compute P1 from P0 + U
            x1 = x0 + ux;
            y1 = y0 + uy;

            % Call the P0P1 version
            idx = obj.getNodeIndicesOnLineP0P1(x0, y0, x1, y1, tol);
        end

        function idx = getNodeIndicesOnRayP0P1(obj, x0, y0, x1, y1, tol)
            % getNodeIndicesOnRayP0P1
            % Returns node indices near the ray starting at P0 = (x0,y0)
            % and going through P1 = (x1,y1).

            if nargin < 6
                tol = obj.gleps;
            end

            p0 = [x0, y0];
            p1 = [x1, y1];

            d = p1 - p0;
            L = norm(d);

            if L < eps
                error('P0 and P1 must be different to define a ray.');
            end

            u = d / L;  % unit vector

            % Relative coordinates of nodes
            rel = obj.nodes(:,1:2) - p0;

            % Projection scalar (along-ray coordinate)
            alpha = rel * u';  % dot product

            % Perpendicular distance from ray axis
            dist = sqrt(sum((rel - alpha*u).^2, 2));

            % Ray condition: alpha >= 0
            mask = alpha >= 0;

            % Distance condition
            idx = find(mask & (dist < tol));
        end

        function idx = getNodeIndicesOnRayP0U(obj, x0, y0, ux, uy, tol)
            % getNodeIndicesOnRayP0U
            % Returns node indices near the ray starting at P0 = (x0,y0)
            % in direction U = (ux,uy).

            % Default tolerance
            if nargin < 6
                tol = obj.gleps;
            end

            % Check direction vector is nonzero
            if ux == 0 && uy == 0
                error('Direction vector U = (ux, uy) must be nonzero.');
            end

            % Compute P1 = P0 + U
            x1 = x0 + ux;
            y1 = y0 + uy;

            % Call the ray version using P0→P1
            idx = obj.getNodeIndicesOnRayP0P1(x0, y0, x1, y1, tol);
        end

        function idx = getNodeIndicesOnSegment(obj, x0, y0, x1, y1, tol)
            % getNodeIndicesOnSegment
            % Returns node indices located on the finite segment from (x0,y0) to (x1,y1)

            % Handle default tolerance
            if nargin < 6
                tol = obj.gleps;
            end

            p0 = [x0, y0];
            p1 = [x1, y1];

            d = p1 - p0;
            L = norm(d);

            if L < eps
                error('Segment endpoints P0 and P1 must be distinct.');
            end

            u = d / L; % Unit direction vector

            % Relative coordinates of nodes from P0
            rel = obj.nodes(:,1:2) - p0;

            % Projection parameter alpha (distance along the line)
            alpha = rel * u';

            % Perpendicular distance from the line
            dist = sqrt(sum((rel - alpha * u).^2, 2));

            % Conditions:
            % 1. Within tolerance of the infinite line
            % 2. Projection is between 0 and the length L (with a small buffer for tol)
            mask = (dist < tol) & (alpha >= -tol) & (alpha <= L + tol);

            idx = find(mask);
        end

        function idx = getNodeIndicesOnEdges(obj,eList)
            switch obj.etype
                case 'QL4'
                    idx = obj.edges(eList,1:4);
                    idx = unique(idx(:));
            end
        end

        function idx = getNodeIndicesInCircle(obj, x0, y0, r, tol)
            % getNodeIndicesInCircle
            % Returns node indices located inside a circle of radius r
            % centered at (x0, y0).
            %
            % A node is considered inside if:
            %     distance(node, centre) <= r + tol

            % Default tolerance
            if nargin < 5
                tol = obj.gleps;
            end

            if r < 0
                error('Radius r must be non‑negative.');
            end

            % Node coordinates
            X = obj.nodes(:,1);
            Y = obj.nodes(:,2);

            % Distance from centre
            dist = hypot(X - x0, Y - y0);

            % Inside test
            idx = find(dist <= r + tol);
        end

        function idx = getNodeIndicesOnCircle(obj, x0, y0, r, tol)
            % getNodeIndicesOnCircle
            % Returns node indices lying on a circle of radius r centered at (x0, y0).
            %
            % A node is considered on the circle if:
            %   |distance(node, centre) - r| < tol

            % Default tolerance
            if nargin < 5
                tol = obj.gleps;
            end

            if r < 0
                error('Radius r must be non‑negative.');
            end

            % Node coordinates
            X = obj.nodes(:,1);
            Y = obj.nodes(:,2);

            % Distance from centre
            dist = hypot(X - x0, Y - y0);

            % Find nodes lying on the circle (within tolerance)
            idx = find(abs(dist - r) < tol);
        end

        function idx = getNodeIndicesOutCircle(obj, x0, y0, r, tol)
            % getNodeIndicesOutCircle
            % Returns node indices located outside a circle of radius r
            % centered at (x0, y0).
            %
            % A node is considered outside if:
            %     distance(node, centre) >= r - tol

            % Default tolerance
            if nargin < 5
                tol = obj.gleps;
            end

            if r < 0
                error('Radius r must be non‑negative.');
            end

            % Node coordinates
            X = obj.nodes(:,1);
            Y = obj.nodes(:,2);

            % Distance from centre
            dist = hypot(X - x0, Y - y0);

            % Outside test
            idx = find(dist >= r - tol);
        end

        function idx = getNodeIndicesOnArcCP0P1(obj, x0, y0, x1, y1, x2, y2, tol)
            % getNodeIndicesOnArcCP0P1
            % Returns node indices on a CCW arc from P1 to P2 centered at (x0,y0)

            % Handle default tolerance
            if nargin < 8
                tol = obj.gleps;
            end

            % Center and reference points
            C = [x0, y0];
            P1 = [x1, y1];
            P2 = [x2, y2];

            % Radius from first point
            R = norm(P1 - C);

            % Get angles of start and end points
            v1 = P1 - C;
            v2 = P2 - C;
            phi1 = atan2(v1(2), v1(1));
            phi2 = atan2(v2(2), v2(1));

            % Normalize phi2 relative to phi1 for a CCW sweep
            % This ensures the arc goes from P1 to P2 in positive direction
            if phi2 < phi1
                phi2 = phi2 + 2*pi;
            end

            % Node coordinates relative to center
            rel = obj.nodes(:, 1:2) - C;

            % 1. Radial distance check
            dist = hypot(rel(:,1), rel(:,2));
            radial_mask = abs(dist - R) < tol;

            % 2. Angular sweep check
            % Get angles of all nodes
            node_phi = atan2(rel(:,2), rel(:,1));

            % We must check node_phi, node_phi + 2*pi, and node_phi - 2*pi
            % to see if any representation falls within [phi1, phi2]
            % Alternatively, shift node_phi to be relative to phi1:
            node_phi_rel = mod(node_phi - phi1, 2*pi);
            sweep_total = phi2 - phi1;

            % Account for numerical precision at boundaries
            angular_mask = (node_phi_rel >= -tol/R) & (node_phi_rel <= sweep_total + tol/R);

            % Combine masks
            idx = find(radial_mask & angular_mask);
        end

        function idx = getNodeIndicesOnArcCP0A(obj, x0, y0, x1, y1, angle, tol)
            % getNodeIndicesOnArcCP0A
            % C = (x0, y0), P0 = (x1, y1), angle = sweep in radians (CCW if positive)

            if nargin < 6, tol = obj.gleps; end

            C = [x0, y0];
            P0 = [x1, y1];

            % Calculate radius and start angle
            v0 = P0 - C;
            R = norm(v0);
            phi_start = atan2(v0(2), v0(1));

            % Node coordinates relative to center
            rel = obj.nodes(:, 1:2) - C;
            dist = hypot(rel(:,1), rel(:,2));

            % Radial mask
            radial_mask = abs(dist - R) < tol;

            % Angular mask
            node_phi = atan2(rel(:,2), rel(:,1));
            % Shift node angles to start at 0 from phi_start
            % Using mod(..., 2*pi) handles the wrap-around
            node_phi_rel = mod(node_phi - phi_start, 2*pi);

            if angle >= 0
                % Counter-clockwise sweep
                angular_mask = (node_phi_rel <= angle + tol/R);
            else
                % Clockwise sweep
                % Map [0, 2pi] to [-2pi, 0]
                node_phi_rel_cw = node_phi_rel - 2*pi;
                angular_mask = (node_phi_rel_cw >= angle - tol/R);
            end

            idx = find(radial_mask & angular_mask);
        end

        function idx = getEdgeIndicesOnLineP0P1(obj, x0, y0, x1, y1, tol)
            % getEdgeIndicesOnLineP0P1
            % Returns indices of edges that lie entirely on the line segment P0->P1.
            % An edge is on the line if both of its nodes are on the line.

            % Handle default tolerance
            if nargin < 6
                tol = obj.gleps;
            end

            % 1. Get the indices of all nodes on this line segment
            node_idx = obj.getNodeIndicesOnLineP0P1(x0, y0, x1, y1, tol);

            % 2. Find edges where BOTH endpoints (node 1 and node 2) are in node_idx.
            % We use a logical mask instead of ismember where possible, or use 'rows'
            % to make ismember fast. For general sets, ismember is correct.
            mask = ismember(obj.edges(:, 1), node_idx) & ismember(obj.edges(:, 2), node_idx);

            % 3. Return the indices of those edges
            idx = find(mask);
        end

        function idx = getEdgeIndicesOnLineP0U(obj, x0, y0, ux, uy, tol)
            % getEdgeIndicesOnLineP0U
            % Returns edge indices whose two end nodes lie on the line through P0
            % in direction U.

            if nargin < 6
                tol = obj.gleps;
            end

            if ux == 0 && uy == 0
                error('Direction vector U = (ux, uy) must be nonzero.');
            end

            nodeIdx = obj.getNodeIndicesOnLineP0U(x0, y0, ux, uy, tol);

            mask = ismember(obj.edges(:,1), nodeIdx) & ...
                ismember(obj.edges(:,2), nodeIdx);

            idx = find(mask);
        end

        function idx = getEdgeIndicesOnCircle(obj, x0, y0, r, tol)
            % getEdgeIndicesOnCircle
            % Returns indices of edges that lie on a circle of radius r
            % centered at (x0, y0).
            %
            % An edge is selected if both of its end nodes are on the circle.

            % Default tolerance
            if nargin < 5
                tol = obj.gleps;
            end

            % 1. Get indices of all nodes lying on this circle boundary
            node_idx = obj.getNodeIndicesOnCircle(x0, y0, r, tol);

            % 2. Find edges where both endpoints are in node_idx
            mask = ismember(obj.edges(:, 1), node_idx) & ...
                ismember(obj.edges(:, 2), node_idx);

            % 3. Return the indices of those edges
            idx = find(mask);
        end

        function idx = getEdgeIndicesInCircle(obj, x0, y0, r, tol)
            % getEdgeIndicesInCircle
            % Returns indices of edges that lie entirely inside (or on)
            % a circle of radius r centered at (x0, y0).
            %
            % An edge is selected if both of its end nodes are within the circle.

            if nargin < 5 || isempty(tol)
                tol = obj.gleps;
            end

            % 1. Extract the start and end node indices for all edges
            edgeNodes = obj.edges(:, 1:2);

            % 2. Retrieve coordinates of all nodes involved in the edges
            p1 = obj.nodes(edgeNodes(:, 1), :);
            p2 = obj.nodes(edgeNodes(:, 2), :);

            % 3. Calculate Euclidean distance from center (x0, y0) to both endpoints
            d1 = hypot(p1(:, 1) - x0, p1(:, 2) - y0);
            d2 = hypot(p2(:, 1) - x0, p2(:, 2) - y0);

            % 4. Select edges where both endpoints are within the radius (with tolerance)
            mask = (d1 <= (r + tol)) & (d2 <= (r + tol));

            % 5. Return the matching edge indices
            idx = find(mask);
        end

        function idx = getEdgeIndicesOutCircle(obj, x0, y0, r, tol)
            % getEdgeIndicesOutCircle
            % Returns indices of edges that lie outside a circle of radius r
            % centered at (x0, y0).
            %
            % Assumes an edge is outside if both endpoints are outside the radius.

            if nargin < 5 || isempty(tol)
                tol = obj.gleps;
            end

            % 1. Extract start and end nodes for all edges
            edgeNodes = obj.edges(:, 1:2);

            % 2. Get coordinates
            p1 = obj.nodes(edgeNodes(:, 1), :);
            p2 = obj.nodes(edgeNodes(:, 2), :);

            % 3. Calculate distance from center to both endpoints
            d1 = hypot(p1(:, 1) - x0, p1(:, 2) - y0);
            d2 = hypot(p2(:, 1) - x0, p2(:, 2) - y0);

            % 4. Select edges where both endpoints are strictly outside the radius
            mask = (d1 >= (r - tol)) & (d2 >= (r - tol));

            % 5. Return indices
            idx = find(mask);
        end

        function idx = getEdgeIndicesOnContact(obj, mz1Name, mz2Name)

            mz1Name = obj.checkMeshZoneExistence(mz1Name);
            mz2Name = obj.checkMeshZoneExistence(mz2Name);

            mask = (ismember(obj.edges(:, 3), obj.mzs.(mz1Name).zi) & ...
                ismember(obj.edges(:, 4), obj.mzs.(mz2Name).zi)) | ...
                (ismember(obj.edges(:, 3), obj.mzs.(mz2Name).zi) & ...
                ismember(obj.edges(:, 4), obj.mzs.(mz1Name).zi));

            idx = find(mask);

        end

    end

end
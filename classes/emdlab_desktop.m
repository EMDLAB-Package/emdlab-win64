classdef emdlab_desktop < handle

    properties
        g = emdlab_g2d_db            % geometry database
        f                           % uifigure
        ax                          % axes
        addPointMode = false        % state flag
    end

    methods

        function obj = emdlab_desktop()

            %% Figure
            obj.f = uifigure( ...
                'Name','emdlab-desktop', ...
                'Position',[100 100 900 600]);

            %% Menus
            mFile = uimenu(obj.f,'Text','File');
            uimenu(mFile,'Text','Open');
            uimenu(mFile,'Text','Save');
            uimenu(mFile,'Text','Exit','Separator','on', ...
                'MenuSelectedFcn',@(s,e)close(obj.f));

            mGeometry = uimenu(obj.f,'Text','Geometry');
            uimenu(mGeometry,'Text','Add Point', ...
                'MenuSelectedFcn',@(s,e)obj.enableAddPoint());

            %% Axes
            obj.ax = uiaxes(obj.f, ...
                'Position',[50 50 800 500]);
            grid(obj.ax,'on')
            axis(obj.ax,'equal')
            hold(obj.ax,'on')
            title(obj.ax,'Geometry Editor')

            %% Mouse click callback
            obj.f.WindowButtonDownFcn = @(s,e)obj.onMouseClick();

        end
    end

    methods (Access = private)

        function enableAddPoint(obj)
            obj.addPointMode = true;
            obj.f.Pointer = 'crosshair';
        end

        function onMouseClick(obj)

            if ~obj.addPointMode
                return
            end

            % Get click position in axes coordinates
            cp = obj.ax.CurrentPoint;
            x = cp(1,1);
            y = cp(1,2);

            % Add point to geometry database
            obj.g.addPoint(x,y);

            % Plot point
            obj.g.showSketch(obj.ax);

            % Reset mode
            obj.addPointMode = false;
            obj.f.Pointer = 'arrow';

        end
    end
end

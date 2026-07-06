function [f,va] = emdlab_r3d_geometry(bgFlag, hflag)

if nargin == 0
    bgFlag = 1;
end

if nargin>2
    hflag = 1;
end

%% ================= FIGURE ==========================================
f = figure( ...
    'Visible', 'off', ...
    'Position', [0 0 800 600], ...
    'Tag', 'main', ...
    'Name', 'EMDLAB', ...
    'Clipping', 'off', ...
    'Color', [0.2 0.2 0.2], ...
    'Renderer', 'opengl', ...
    'Units', 'pixels', ...
    'ToolBar','none', ...
    'DockControls','off',...
    'MenuBar','none');

movegui(f, 'center');

f.UserData = struct( ...
    'cp', [0 0], ...
    'cp0', [0 0], ...
    'isMouseDown', false, ...
    'isMouseLeft', false, ...
    'isMouseRight', false, ...
    'isDragging', false, ...
    'zf', 1, ...
    'hoverPatch', gobjects(0), ...
    'dragThreshold', 3);

%% ================= MAIN VIEW AXIS ==================================
va = axes(f, ...
    'Tag', 'va', ...
    'Units', 'normalized', ...
    'Position', [0 0 1 1], ...
    'Clipping', 'off');

va.Toolbar.Visible = 'off';

camzoom(va,0.7);
axis(va,'off','vis3d');

%% ================= BACKGROUND ======================================
bg = [];
if bgFlag
    bg = axes(f, ...
        'Tag', 'bg', ...
        'Units', 'normalized', ...
        'Position', [0 0 1 1], ...
        'HandleVisibility', 'off', ...
        'Clipping', 'off', ...
        'HitTest', 'off', ...
        'PickableParts', 'none');
    bg.Toolbar.Visible = 'off';

    uistack(bg,'bottom');
    setBackground(bg);
    %     drawBackgroundLabel(bg, 'emdlab_win64');
end

%% ================= SMALL COORDINATE AXIS ===========================
ca = axes(f, ...
    'Tag', 'ca', ...
    'Units', 'pixels', ...
    'Position', [20 20 100 100], ...
    'NextPlot','add');

ca.Toolbar.Visible = 'off';
if ~verLessThan('matlab','9.5')
    ca.Toolbar = [];
end
ca.Interactions = [];
ca.PickableParts = 'none';

axis(ca,'off','vis3d');

drawAxisArrow(ca, [1 0 0], 'r', 'X');
drawAxisArrow(ca, [0 1 0], 'b', 'Y');
drawAxisArrow(ca, [0 0 1], [0 0.6 0], 'Z');

patch(ca, ...
    'Faces', [1 2; 3 4; 5 6], ...
    'Vertices', [-1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1], ...
    'EdgeColor','none');

ca.HandleVisibility = 'off';

if bgFlag
    f.UserData.txt = drawBackgroundLabel(bg,'');
else
    f.UserData.txt = '';
end

va.UserData.camMovH = [0,0,0];
va.UserData.camPosH = [0,0,0];
va.UserData.rFlag = true;

%% ================= CALLBACKS =======================================
f.WindowButtonDownFcn   = @fButtonDownFcn;
f.WindowButtonMotionFcn = @fButtonMotionFcn;
f.WindowButtonUpFcn     = @fButtonUpFcn;
f.WindowScrollWheelFcn  = @fScrollFcn;
f.WindowKeyPressFcn     = @fKeyFcn;
f.SizeChangedFcn        = @fResizeFcn;

guidata(f, guihandles(f));
f.Visible = 'on';

%% ================= CALLBACK FUNCTIONS ===============================

    function fButtonDownFcn(~,~)
        ud = f.UserData;
        ud.isMouseDown = true;
        ud.cp = f.CurrentPoint;
        ud.cp0 = ud.cp;
        ud.isMouseLeft = strcmpi(f.SelectionType,'normal');
        ud.isMouseRight = strcmpi(f.SelectionType,'alt');
        ud.isDragging = false;
        f.UserData = ud;
    end

    function fButtonMotionFcn(src,~)
        ud = f.UserData;

        % ---------- DRAG MODE: rotate only ----------
        if ud.isMouseDown && (ud.isMouseLeft || ud.isMouseRight)
            cpNew = f.CurrentPoint;
            dxy = cpNew - ud.cp;

            if norm(cpNew - ud.cp0) > ud.dragThreshold
                ud.isDragging = true;
            end

            ud.cp = cpNew;

            if ud.isDragging
                dx = dxy(1);
                dy = dxy(2);

                if ud.isMouseLeft
                    rotScale = 0.35;
                    camorbit(va, -rotScale*dx, -rotScale*dy, 'camera');
                    camorbit(ca, -rotScale*dx, -rotScale*dy, 'camera');
                else
                    mousePanVa(dx,dy);
                end
            end

            f.UserData = ud;
            drawnow limitrate;
            return
        end

        % ---------- HOVER MODE: highlight only ----------
        hoverHighlight(src);
    end

    function fButtonUpFcn(src,~)
        ud = f.UserData;

        % click selection only if it was not a drag
        if ud.isMouseLeft && ~ud.isDragging
            h = hittest(src);
            if isa(h,'matlab.graphics.primitive.Patch') && ancestor(h,'axes') == va
                e.Button = 1;
                emdlab_flib_selectPatchCallbackGM(h,e);
            end
        end

        ud.isMouseDown = false;
        ud.isMouseLeft = false;
        ud.isMouseRight = false;
        ud.isDragging = false;
        f.UserData = ud;
    end

    function fScrollFcn(~,eventData)
        ud = f.UserData;

        if eventData.VerticalScrollCount < 0
            camzoom(va,0.9);
            ud.zf = ud.zf * 0.9;
        else
            camzoom(va,1.1);
            ud.zf = ud.zf * 1.1;
        end

        f.UserData = ud;
    end

    function fKeyFcn(~,eventData)
        key = lower(eventData.Key);
        mods = lower(string(eventData.Modifier));

        hasNoMod = isempty(mods);
        hasShift = numel(mods)==1 && mods=="shift";
        hasCtrl  = numel(mods)==1 && mods=="control";

        switch key
            %% ---------- ORBIT (no modifier) ----------
            case 'rightarrow'
                if hasNoMod
                    camorbit(ca,-10,0,'camera');
                    camorbit(va,-10,0,'camera');
                elseif hasShift
                    camroll(va,-5);
                    camroll(ca,-5);
                elseif hasCtrl
                    keyboardPanVa(+10,0);
                end

            case 'leftarrow'
                if hasNoMod
                    camorbit(ca,10,0,'camera');
                    camorbit(va,10,0,'camera');
                elseif hasShift
                    camroll(va,5);
                    camroll(ca,5);
                elseif hasCtrl
                    keyboardPanVa(-10,0);
                end

            case 'uparrow'
                if hasNoMod
                    camorbit(ca,0,-10,'camera');
                    camorbit(va,0,-10,'camera');
                elseif hasShift
                    ud = f.UserData;
                    camzoom(va,1.1);
                    ud.zf = ud.zf*1.1;
                    f.UserData = ud;
                elseif hasCtrl
                    keyboardPanVa(0,+10);
                end

            case 'downarrow'
                if hasNoMod
                    camorbit(ca,0,10,'camera');
                    camorbit(va,0,10,'camera');
                elseif hasShift
                    ud = f.UserData;
                    camzoom(va,0.9);
                    ud.zf = ud.zf*0.9;
                    f.UserData = ud;
                elseif hasCtrl
                    keyboardPanVa(0,-10);
                end

                %% ---------- PRESET VIEWS ----------
            case 'z'
                if hasNoMod
                    view(va,[0 0 1]);
                    view(ca,[0 0 1]);
                end
            case 'x'
                if hasNoMod
                    view(va,[1 0 0]);
                    view(ca,[1 0 0]);
                end
            case 'y'
                if hasNoMod
                    view(va,[0 1 0]);
                    view(ca,[0 1 0]);
                end
            case 'i'
                if hasNoMod
                    camorbit(ca,30,20,'camera');
                    camorbit(va,30,20,'camera');
                end

                %% ---------- ROLL ----------
            case 'space'
                if hasNoMod
                    camroll(va,90);
                    camroll(ca,90);
                end

                %% ---------- SAVE EPS ----------
            case 's'
                if hasCtrl
                    if ~isempty(bg), bg.Visible='off'; end
                    f.Renderer='painters';
                    saveas(f,'view.eps');
                    if ~isempty(bg), bg.Visible='on'; end
                end

                %% ---------- ZOOM FIT ----------
            case 'f'
                if hasNoMod
                    ud = f.UserData;
                    camzoom(va,1/ud.zf);
                    ud.zf = 1;
                    f.UserData = ud;
                end
            case 'r'

                va.CameraPosition = va.UserData.camMovH;
                va.CameraTarget   =  va.UserData.camPosH;
                va.UserData.rFlag = true;

        end
    end

    function keyboardPanVa(dx,dy)
        campos_ = va.CameraPosition;
        camtar_ = va.CameraTarget;
        camup_  = va.CameraUpVector;

        forward = (camtar_ - campos_);
        forward = forward ./ norm(forward);

        right = cross(forward, camup_);
        right = right ./ norm(right);

        up = camup_ ./ norm(camup_);

        panSpeed = 0.2;
        move = (dx * right + dy * up) * panSpeed;

        va.CameraPosition = campos_ - move;
        va.CameraTarget   = camtar_ - move;

        % ca intentionally NOT panned
    end


    function mousePanVa(dx,dy)


        campos_ = va.CameraPosition;
        camtar_ = va.CameraTarget;
        camup_  = va.CameraUpVector;

        % Camera basis vectors
        forward = camtar_ - campos_;
        dist = norm(forward);
        forward = forward/dist;

        right = cross(forward,camup_);
        right = right/norm(right);

        up = cross(right,forward);
        up = up/norm(up);

        % Translation proportional to camera distance
        s = dist/9800;

        move = (-dx*right - dy*up)*s;



        if va.UserData.rFlag
            va.UserData.camMovH = campos_;
            va.UserData.camPosH = camtar_;
            va.UserData.rFlag = false;
        end

        va.CameraPosition = campos_ + move;
        va.CameraTarget   = camtar_ + move;

    end
    function fResizeFcn(~,~)
        handles = guihandles(f);
        if isfield(handles,'cba')
            handles.cba.Units='normalized';
            handles.cba.Position(2)=0.2;
            handles.cba.Position(4)=0.7;
            handles.cba.Units='pixels';
            handles.cba.Position(1)=30;
            handles.cba.Position(3)=20;
        end
    end

%% ================= HELPERS =========================================

    function hoverHighlight(src)
        if ~hflag; return; end
        ud = f.UserData;
        h = hittest(src);

        % Restore previously hovered patch
        if ~isempty(ud.hoverPatch) && isgraphics(ud.hoverPatch)
            if isa(ud.hoverPatch,'matlab.graphics.primitive.Patch')
                restorePatchAppearance(ud.hoverPatch);
            end
        end

        ud.hoverPatch = gobjects(0);

        % Default title when not hovering a patch
        f.UserData.txt.String = '';

        % Highlight current hovered patch only if it belongs to va
        if isa(h,'matlab.graphics.primitive.Patch') && ancestor(h,'axes') == va
            set(h, 'FaceColor', 'y', 'FaceAlpha', 1, 'edgecolor', [0.2,0.2,0.2]);
            ud.hoverPatch = h;

            f.UserData.txt.String = h.UserData.Tag;

        end

        f.UserData = ud;
        drawnow limitrate;
    end


    function restorePatchAppearance(hp)
        if ~isgraphics(hp)
            return
        end

        % default restore state
        set(hp, 'FaceColor', hp.UserData.c, 'FaceAlpha', 1, 'edgecolor', [0.2,0.2,0.2]);
    end

    function setBackground(ax)
        n = 250;
        d = n/3;

        c1 = [1 1 1; d^2 d 1; n^2 n 1]\[1;0.80;1];
        c2 = [1 1 1; d^2 d 1; n^2 n 1]\[1;0.81;1];
        c3 = [1 1 1; d^2 d 1; n^2 n 1]\[1;0.82;1];

        r = (1:n)';
        m = zeros(n,n,3);

        m(:,:,1) = repmat(c1(1)*r.^2 + c1(2)*r + c1(3),1,n);
        m(:,:,2) = repmat(c2(1)*r.^2 + c2(2)*r + c2(3),1,n);
        m(:,:,3) = repmat(c3(1)*r.^2 + c3(2)*r + c3(3),1,n);

        imagesc(ax,m);
        axis(ax,'off');
    end

    function txt = drawBackgroundLabel(ax,str)
        txt=text(ax,0.5,0.95,str, ...
            'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontName','Consolas', ...
            'FontSize',24, ...
            'FontWeight','bold', ...
            'Color',[1 1 1]*0.2, ...
            'Interpreter','none', ...
            'HitTest','off');
    end

    function drawAxisArrow(ax,dir,color,labelText)
        emdlab_flib_arrow3d([0 dir(1)], [0 dir(2)], [0 dir(3)], ...
            0.85,0.045,0.1,color,ax);
        text(ax,1.1*dir(1),1.1*dir(2),1.1*dir(3), ...
            labelText,'Color',color,'FontSize',10);
    end

end

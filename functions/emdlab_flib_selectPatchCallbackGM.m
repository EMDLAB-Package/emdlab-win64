function emdlab_flib_selectPatchCallbackGM(h, e)

    if e.Button == 1

        if isequal('FaceAlpha', 1)
            return;
        end

        % view axis
        ax = h.Parent;

        % patch
        plt = ax.findobj('FaceColor', 'r');

        if ~isempty(plt)
            set(plt, 'FaceColor', h.UserData.color);
        end

        set(h, 'FaceColor', 'r');
        title(ax, h.UserData.title,'Interpreter','none');
        set(gcf, 'Name', h.UserData.title)
        
        drawnow;
    end

end

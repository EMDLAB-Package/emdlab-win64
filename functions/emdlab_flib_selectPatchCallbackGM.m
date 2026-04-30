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
            set(plt, 'FaceColor', 'c');
        end

        set(h, 'FaceColor', 'r');
        title(ax, h.UserData,'Interpreter','none');
        set(gcf, 'Name', h.UserData)
        drawnow;
    end

end

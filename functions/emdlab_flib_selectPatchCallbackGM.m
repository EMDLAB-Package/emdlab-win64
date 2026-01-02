function emdlab_flib_selectPatchCallbackGM(h, e)

    if e.Button == 1

        if isequal('FaceAlpha', 1)
            return;
        end

        % view axis
        ax = h.Parent;

        % patch
        plt = ax.findobj('FaceAlpha', 1);

        if ~isempty(plt)
            set(plt, 'FaceColor', 'c', 'FaceAlpha', 0.8);
        end

        set(h, 'FaceColor', 'r', 'FaceAlpha', 1);
        title(ax, h.UserData);
        set(gcf, 'Name', h.UserData)
        drawnow;
    end

end

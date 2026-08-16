function emdlab_flib_selectPatchCallbackGM1(h, e)

    if e.Button == 1

        % view axis
        ax = h.Parent;

        % patch
        for i = 1:numel(ax.Children)
            plt = ax.Children(i);
            if ~isempty(plt.UserData)
            if plt.UserData.wasActive   
                if ~isequal(h,plt)
                plt.UserData.wasActive = false;
                set(plt, 'FaceColor', plt.UserData.color);
                break;
                end
            end
            end
        end

        set(h, 'FaceColor', 'y');
        title(ax, h.UserData.title,'Interpreter','none');
        set(gcf, 'Name', h.UserData.title)
        
        drawnow;
    end

end

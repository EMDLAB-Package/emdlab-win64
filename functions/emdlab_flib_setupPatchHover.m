function emdlab_flib_setupPatchHover(fig)

    patches = findobj(fig, 'Type', 'patch');

    for k = 1:numel(patches)
        patches(k).UserData.originalFaceColor = patches(k).FaceColor;
    end

    previousPatch = gobjects(1);

    % Store original title for each axes
    axesList = findobj(fig, 'Type', 'axes');
    originalTitles = cell(size(axesList));

    for k = 1:numel(axesList)
        originalTitles{k} = axesList(k).Title.String;
    end

    fig.WindowButtonMotionFcn = @hoverFcn;


    function hoverFcn(src, ~)

        obj = hittest(src);

        if isgraphics(obj) && strcmp(obj.Type, 'patch')

            currentPatch = obj;

            if ~isequal(currentPatch, previousPatch)

                % Restore previous patch
                if isgraphics(previousPatch)
                    previousPatch.FaceColor = ...
                        previousPatch.UserData.originalFaceColor;
                end

                % Highlight current patch
                currentPatch.FaceColor = 'yellow';

                previousPatch = currentPatch;

                % Show patch title on axes
                ax = currentPatch.Parent;

                title(ax, currentPatch.UserData.title, ...
                    'Interpreter', 'none');

            end

        else

            % Restore previous patch
            if isgraphics(previousPatch)

                ax = previousPatch.Parent;

                previousPatch.FaceColor = ...
                    previousPatch.UserData.originalFaceColor;

                % Restore original axes title
                idx = find(axesList == ax, 1);

                if ~isempty(idx)
                    title(ax, originalTitles{idx}, ...
                        'Interpreter', 'none');
                end

                previousPatch = gobjects(1);
            end

        end
    end
end
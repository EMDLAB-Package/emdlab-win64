classdef emdlab_ui_console < handle

    properties

        % flag to print the elapsed times
        printFlag (1,1) logical = true;

    end

    methods

        function setPrintFlag(obj, newValue)
            obj.printFlag = newValue;
        end

        function timeHolder = dispLine(obj)
            if obj.printFlag
                disp('-------------------------------------------------------');
            end
            timeHolder = tic;
        end

        function dispMessage(obj, txt, timeHolder)
            if obj.printFlag
                disp(txt);
                toc(timeHolder);
            end
        end

        function dispMessageLine(obj, txt, timeHolder)
            if obj.printFlag
                disp(txt);
                toc(timeHolder);
                disp('-------------------------------------------------------');
            end
        end

        function fprintf(obj, varargin)
            if obj.printFlag
                fprintf(varargin{:});
            end
        end

    end

end
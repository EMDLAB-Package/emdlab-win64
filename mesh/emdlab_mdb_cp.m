% EMDLAB: Electrical Machines Design Laboratory
% common properties for all mesh database classes

classdef emdlab_mdb_cp < handle

    properties

        % mesh zones
        mzs (1,1) struct;

        % materials
        mts (1,1) struct;

        % states
        isGlobalMeshGenerated (1,1) logical = false;

        % moving regions
        mrgs (1,1) struct;

        % interface mesh zones: interface between static and moving regions
        imzs (1,1) struct;

        % a vector to store material names
        materialNames (1,:) string;

    end

    properties (Dependent = true)

        % number of mesh zones
        Nmzs (1,1) double {mustBeNonnegative, mustBeInteger};

        % number of materials
        Nmts (1,1) double {mustBeNonnegative, mustBeInteger};

    end

    methods

        function obj = emdlab_mdb_cp()
        end

        function y = get.Nmzs(obj)
            y = numel(fieldnames(obj.mzs));
        end

        function y = get.Nmts(obj)
            y = numel(fieldnames(obj.mts));
        end

        function y = getMeshZoneNames(obj)
            y = string(fieldnames(obj.mzs))';
        end

        function y = getDefaultMeshZoneName(obj)

            index = 0;
            mzNames = fieldnames(obj.mzs);

            while true
                index = index + 1;
                y = ['Zone', num2str(index)];
                if ~ismember(y, mzNames), break; end
            end

        end

        function mzName = checkMeshZoneExistence(obj, mzName)

            mzName = erase(mzName, ' ');
            if ~isfield(obj.mzs, mzName)
                error('Specified mesh zone does not exist.');
            end

        end

        function mzName = checkMeshZoneNonExistence(obj, mzName)

            mzName = erase(mzName, ' ');
            if isfield(obj.mzs, mzName)
                error('Specified mesh zone already exist.');
            end

        end

        function setMeshZoneColor(obj, mzName, R, G, B)

            if ischar(mzName)
                mzName = obj.checkMeshZoneExistence(mzName);
                obj.mzs.(mzName).color = [R,G,B]/255;
            elseif isstring(mzName)
                for i = 1:numel(mzName)
                    mzNameChar = obj.checkMeshZoneExistence(mzName(i));
                    obj.mzs.(mzNameChar).color = [R,G,B]/255;
                end
            else
                error('Wrong input type, mzName must be char or string');
            end

        end

        function setmzc(varargin)

            setMeshZoneColor(varargin{:});

        end

        function addMeshZone(obj, varargin)

            if nargin == 2

%                 if ~isa(varargin{1}, 'emdlab_m2d_tmz')
%                     error('Mesh zone class must be <emdlab_m2d_tmz>.');
%                 end

                mzName = obj.getDefaultMeshZoneName;
                mzptr = varargin{1};

            elseif nargin == 3
                mzName = obj.checkMeshZoneNonExistence(varargin{1});

%                 if ~isa(varargin{2}, 'emdlab_m2d_tmz')
%                     error('Mesh zone class must be <emdlab_m2d_tmz>.');
%                 end

                mzptr = varargin{2};
            else
                error('Wrong number of arguments.');
            end

            % adding new mesh zone
            obj.mzs.(mzName) = mzptr;
            obj.mzs.(mzName).material = 'air';
            obj.mzs.(mzName).color = rand(1,3);

            % changing states
            obj.makeFalse_isGlobalMeshGenerated;

        end

        function addmz(obj, varargin)

            obj.addMeshZone(varargin{:});

        end

        function removeMeshZone(obj, mzName)

            mzName = obj.checkMeshZoneExistence(mzName);
            % delete(obj.mzs.(mzName));
            obj.mzs = rmfield(obj.mzs, mzName);
            % changing states
            obj.makeFalse_isGlobalMeshGenerated;

        end

        function removemz(obj, varargin)

            obj.removeMeshZone(varargin{:});

        end

        function changeMeshZoneName(obj, mzName, newName)

            newName = obj.checkMeshZoneNonExistence(newName);
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(newName) = copy(obj.mzs.(mzName));
            obj.mzs = rmfield(obj.mzs, mzName);

        end

        function clearAllmzs(obj)

            obj.mzs = struct();

        end

        % moving region
        function mrName = checkMovingRegionExistence(obj, mrName)

            mrName = erase(mrName, ' ');

            if ~isfield(obj.mrgs, mrName)
                error('Specified moving region does not exist.');
            end

        end

        function mrName = checkMovingRegionNonExistence(obj, mrName)

            mrName = erase(mrName, ' ');

            if isfield(obj.mzs, mrName)
                error('Specified moving region name already exist.');
            end

        end

        function defineMovingRegion(obj, mrName, varargin)

            mrName = obj.checkMovingRegionNonExistence(mrName);
            obj.mrgs.(mrName) = emdlab_solvers_mt2d_star();
            obj.starConnections.(starConnectionName).sci = obj.NstarConnections;

            if all(cellfun(@ischar, varargin))
                coilNames = string(varargin);
            else
                coilNames = string([varargin{:}]);
            end

        end

        function moveRegion(obj, mrName, shiftX, shiftY)
        end

        function rotateRegion(obj, mrName, rotAngle, xc, yc)
        end

        function showMovingRegion(obj, mrName)
        end

        % moving region
        function mzName = checkInterfaceExistence(obj, mzName)

            mzName = erase(mzName, ' ');

            if ~isfield(obj.imzs, mzName)
                error('This interface does not exist.');
            end

        end

        function mzName = checkInterfaceNonExistence(obj, mzName)

            mzName = erase(mzName, ' ');

            if isfield(obj.mzs, mzName)
                error('Specified mesh zone is already defined as interface.');
            end

        end

        function addInterface(obj, mzName, interfaceObject)

            mzName = obj.checkMeshZoneExistence(mzName);
            obj.imzs.(mzName) = interfaceObject;

        end

        function showInterface(obj, mrName)
        end

        % material library
        function setMaterial(obj, mzName, mName)

            mzName = obj.checkMeshZoneExistence(mzName);
            mName = erase(mName, ' ');

            if ~ismember(mName, fieldnames(obj.mts))
                error(['Material <<', mName, '>> does not found.']);
            end

            obj.mzs.(mzName).material = mName;

        end

        function materialName = getMaterialByZoneIndex(obj, zi)

            mzNames = fieldnames(obj.mzs);
            for i = 1:obj.Nmzs
                if obj.mzs.(mzNames{i}).zi == zi
                    materialName = obj.mzs.(mzNames{i}).material;
                    return;
                end
            end
            materialName = '';

        end

        function addMaterial(obj, matrialName, mObject)

            if nargin == 2
                obj.mts.(matrialName) = eval(['emdlab_mlib_', matrialName]);
            else
                obj.mts.(matrialName) = mObject;
            end

        end

    end

end
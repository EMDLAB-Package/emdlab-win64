% EMDLAB: Electrical Machines Design Laboratory
% two-dimensional magnetic-transient motion definition
% motion can be linear or rotational

classdef emdlab_solvers_mt2d_motion < handle & matlab.mixin.SetGet

    properties

        % coil arm names: all coil arms are solid
        meshZones (1,:) string;

        % motion index
        mi (1,1) double;

        % interface mesh zone
        interface (1,:) char;

        % motion history -> [xsh, ysh, x0, y0, theta]
        motionHistory (:,5) double;

        idx

    end

end
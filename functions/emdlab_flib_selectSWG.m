% **********************************************************************
% this function gets the desired area for conductor and the number of
% strands, and it return a row vector of wire diameters and corresponding
% gauges, swd: select wire diameters
% developer: Ali Jamali-Fard, ComProgExpert.com
% **********************************************************************
function wireD = emdlab_flib_selectSWG(desiredArea, Nstrands)

% available wire diameters (swg)
d = [12.7,11.786,10.973,10.16,9.449,8.839,8.23,7.62,7.01,6.401,5.893,5.385,4.877,4.47,4.064 ...
,3.658,3.251,2.946,2.642,2.337,2.032,1.829,1.626,1.422,1.219,1.016,0.914,0.813,0.711,0.61 ...
,0.559,0.508,0.4572,0.4166,0.3759,0.3454,0.315,0.2946,0.2743,0.254,0.2337,0.2134,0.193,0.1727 ...
,0.1524,0.1321,0.1219,0.1118,0.1016,0.0914,0.0813,0.0711,0.061,0.0508,0.0406,0.0305,0.0254];

% sort wire diameters
d = sort(d);

% area of wires (round wires)
a = pi*d.^2/4;

for i = 1:length(d)
    if a(i)*Nstrands>desiredArea
        break;
    end
end
wireD = d(i)*ones(1,Nstrands);

end
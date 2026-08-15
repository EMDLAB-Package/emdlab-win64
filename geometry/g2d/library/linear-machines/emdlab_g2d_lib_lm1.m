% linear machine models

function m = emdlab_g2d_lib_lm1(g, p, q, Nwl, shp, ws, wt, ds, dymc, Lemc, mcta, dar, dyrc, gap, np, Ley)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    q = 2;
    p = 4;
    Nwl = 2;
    shp = 1;
    ws = 10;
    wt = 8;
    ds = 40;
    dymc = 18;
    dar = 4;
    dyrc = 8;
    gap = 10;
    np = 2;
    Lemc = 30;
    mcta = 15;
    Ley = 30;
end

L2 = Lemc + (ds + dymc) * tan(deg2rad(mcta));
tau_s = ws + wt;

if Nwl == 1

    % number of slots in moving core
    Ns = 3*q*p;

    x_tmp = [wt/2;wt/2;wt/2+ws;wt/2+ws];
    y_tmp = gap/2 + [0;ds;ds;0];
    tmp = 0:tau_s:(Ns-1)*tau_s;
    tmp = repmat(tmp,4,1);
    x_tmp = repmat(x_tmp,1,Ns);
    x_tmp = x_tmp + tmp;
    y_tmp = repmat(y_tmp,1,Ns);

    x_tmp = x_tmp(:)';
    y_tmp = y_tmp(:)';

    x_tmp = [0,x_tmp,Ns*tau_s,Ns*tau_s,0];
    y_tmp = [gap/2,y_tmp,gap/2,gap/2+ds+dymc,gap/2+ds+dymc];

    x0 = Ns*tau_s/2;
    x_tmp = x_tmp - x0;
    
    x_tmp(end) = x_tmp(end) - L2;
    x_tmp(end-1) = x_tmp(end-1) + L2;
    x_tmp(1) = x_tmp(1) - Lemc;
    x_tmp(end-2) = x_tmp(end-2) + Lemc;

    l_mc = g.addClosedPolylineLoop(x_tmp,y_tmp);

    for i = 1:Ns
        idx1 = 2 + (i-1)*4;
        eidx_new1 = g.addSegmentByCoordinates(x_tmp(idx1),y_tmp(idx1),x_tmp(idx1+3),y_tmp(idx1+3));
        g.addFace('ca_' + string(i), g.addLoop(g.getEdgeLeftLoop(eidx_new1)));
        g.setFaceColor('ca_' + string(i),0,255,255);
    end

    x1 = np*Ns*tau_s/2;
    y1 = gap/2+ds+dymc+Ley;
    ei = g.addSegmentByCoordinates(x0+Lemc,gap/2,x1,gap/2);
    g.addSegmentByCoordinates(x1,gap/2,x1,y1);
    g.addSegmentByCoordinates(x1,y1,-x1,y1);
    g.addSegmentByCoordinates(-x1,y1,-x1,gap/2);
    g.addSegmentByCoordinates(-x1,gap/2,-x0-Lemc,gap/2);
    g.addFace('moving_air', g.addLoop(g.getEdgeLeftLoop(ei)))
    g.setFaceColor('moving_air', 0, 255, 255)

elseif Nwl == 2

     % number of slots in moving core
    Ns = 3*q*p + 3*q - shp;

    x_tmp = [wt/2;wt/2;wt/2;wt/2+ws;wt/2+ws;wt/2+ws];
    y_tmp = gap/2 + [0;ds/2;ds;ds;ds/2;0];
    tmp = 0:tau_s:(Ns-1)*tau_s;
    tmp = repmat(tmp,6,1);
    x_tmp = repmat(x_tmp,1,Ns);
    x_tmp = x_tmp + tmp;
    y_tmp = repmat(y_tmp,1,Ns);

    x_tmp = x_tmp(:)';
    y_tmp = y_tmp(:)';

    x_tmp = [0,x_tmp,Ns*tau_s,Ns*tau_s,0];
    y_tmp = [gap/2,y_tmp,gap/2,gap/2+ds+dymc,gap/2+ds+dymc];

    x0 = Ns*tau_s/2;
    x_tmp = x_tmp - x0;
    
    x_tmp(end) = x_tmp(end) - L2;
    x_tmp(end-1) = x_tmp(end-1) + L2;
    x_tmp(1) = x_tmp(1) - Lemc;
    x_tmp(end-2) = x_tmp(end-2) + Lemc;

    l_mc = g.addClosedPolylineLoop(x_tmp,y_tmp);

    for i = 1:Ns
        idx1 = 2 + (i-1)*6;
        eidx_new1 = g.addSegmentByCoordinates(x_tmp(idx1),y_tmp(idx1),x_tmp(idx1+5),y_tmp(idx1+5));
        eidx_new2 = g.addSegmentByCoordinates(x_tmp(idx1+1),y_tmp(idx1+1),x_tmp(idx1+4),y_tmp(idx1+4));
        g.addFace('ca_1_' + string(i), g.getEdgeLeftLoop(eidx_new1));
        g.setFaceColor('ca_1_' + string(i),0,255,255);
        g.addFace('ca_2_' + string(i), g.getEdgeLeftLoop(eidx_new2));
        g.setFaceColor('ca_2_' + string(i),0,255,255);
    end

    x1 = np*Ns*tau_s/2;
    y1 = gap/2+ds+dymc+Ley;
    ei = g.addSegmentByCoordinates(x0+Lemc,gap/2,x1,gap/2);
    g.addSegmentByCoordinates(x1,gap/2,x1,y1);
    g.addSegmentByCoordinates(x1,y1,-x1,y1);
    g.addSegmentByCoordinates(-x1,y1,-x1,gap/2);
    g.addSegmentByCoordinates(-x1,gap/2,-x0-Lemc,gap/2);
    g.addFace('moving_air', g.getEdgeLeftLoop(ei))
    g.setFaceColor('moving_air', 0, 255, 255)

else
    error('Number of winding layers must be 1 or 2.');    
end

% add rail parts
ei1 = g.addSegmentByCoordinates(-x1,-gap/2-dyrc-dar-Ley,x1,-gap/2-dyrc-dar-Ley);
ei2 = g.addSegmentByCoordinates(-x1,-gap/2-dyrc-dar,x1,-gap/2-dyrc-dar);
ei3 = g.addSegmentByCoordinates(-x1,-gap/2-dar,x1,-gap/2-dar);
g.addSegmentByCoordinates(-x1,-gap/2,x1,-gap/2);

g.addSegmentByCoordinates(-x1,-gap/2-dyrc-dar-Ley,-x1,-gap/2-dyrc-dar);
g.addSegmentByCoordinates(-x1,-gap/2-dyrc-dar,-x1,-gap/2-dar);
g.addSegmentByCoordinates(-x1,-gap/2-dar,-x1,-gap/2);

g.addSegmentByCoordinates(x1,-gap/2-dyrc-dar-Ley,x1,-gap/2-dyrc-dar);
g.addSegmentByCoordinates(x1,-gap/2-dyrc-dar,x1,-gap/2-dar);
g.addSegmentByCoordinates(x1,-gap/2-dar,x1,-gap/2);

% add faces
g.addFace('moving_core', l_mc);
g.addFace('rail_air', g.getEdgeLeftLoop(ei1));
g.addFace('rail_core', g.getEdgeLeftLoop(ei2));
g.addFace('rail_caluminum', g.getEdgeLeftLoop(ei3));

% set face colors
g.setFaceColor('moving_core',200,200,200)
g.setFaceColor('rail_core',200,200,200)
g.setFaceColor('rail_caluminum',37,177,76);
g.setFaceColor('rail_air',0,255,255)

% visualizations for debug
close all;
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

g.setMeshMaxLength(ds/3);
m = g.generateMesh('mg0');
end
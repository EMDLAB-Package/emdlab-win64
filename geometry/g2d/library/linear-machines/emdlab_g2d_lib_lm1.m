% linear machine models

function emdlab_g2d_lib_lm1(g, parameters)

% defult arguments for debug
if nargin == 0
    g = emdlab_g2d_db;
    q = 2;
    p = 4;
    Nwl = 1;
    ws = 10;
    wt = 8;
    ds = 40;
    dymc = 20;
    gap = 10;
    np = 3;
end

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

    x_tmp = x_tmp + 100;

    L1 = 50;
    L2 = 70;
    x_tmp(end) = x_tmp(end) - L2;
    x_tmp(end-1) = x_tmp(end-1) + L2;
    x_tmp(1) = x_tmp(1) - L1;
    x_tmp(end-2) = x_tmp(end-2) + L1;

    l_mc = g.addClosedPolylineLoop(x_tmp,y_tmp);

    eidx = [4,3,2];
    elist = zeros(1,2*Ns);
    for i = 1:Ns
        idx1 = 2 + (i-1)*4;
        eidx_new = g.addSegmentByCoordinates(x_tmp(idx1),y_tmp(idx1),x_tmp(idx1+3),y_tmp(idx1+3));
        g.addFace('ca_' + string(i), g.addLoop(eidx_new, -(eidx + 4*(i-1))));
        g.setFaceColor('ca_' + string(i),255,127,39);
        elist(2*i-1) = eidx_new;
        elist(2*i) = 1 + 4*i;
    end

    La = Ns*tau_s;
    l2 = g.addRectangleLoop(0,0,np*La,gap+ds+dymc+5);
    elist = [1,elist,4*Ns+2,4*Ns+3,4*Ns+4];
    g.addFace('region_1', l2, g.addLoop(elist))


else

end

% % rotor pole pitch angle
% alpha_p = 2*pi/p;
% 
% % limit embrace
% embrace = min(embrace,0.95);
% embrace = max(embrace,0.05);
% 
% % magnet arc angle
% alpha_m = embrace * alpha_p;
% 
% x1 = Dsh/2;
% y1 = 0;
% 
% x2 = ISD/2 - dm - gap;
% y2 = 0;
% 
% x3 = ISD/2 - gap;
% y3 = 0;
% 
% x4 = x3 + gap/Nagl;
% y4 = 0;
% 
% x5 = x2;
% y5 = x2 * tan(alpha_m/2);
% 
% y5_max = sqrt((ISD/2-gap)^2 - x2^2); 
% y5_flag = false;
% if y5 > y5_max, y5_flag = true; y5= y5_max; end
% 
% [x6,y6] = g.getIntersectionRayCircle(x5,y5,1,0,0,0,ISD/2-gap*gg);
% 
% x7 = (x4) * cos(alpha_p/2);
% y7 = (x4) * sin(alpha_p/2);
% 
% r_tmp = norm([x5,y5]);
% x8 = r_tmp * cos(alpha_p/2);
% y8 = r_tmp * sin(alpha_p/2);
% 
% x9 = (Dsh/2) * cos(alpha_p/2);
% y9 = (Dsh/2) * sin(alpha_p/2);
% 
% if isempty(x6) || y5_flag
% 
%     xc = (x5^2+y5^2-x3^2)/(2*(x5-x3));    
%     p1 = g.addPoint(x1,y1);
%     p2 = g.addPoint(x2,y2);
%     p3 = g.addPoint(x3,y3);
%     p4 = g.addPoint(x4,y4);
%     p5 = g.addPoint(x5,y5);
%     p6 = g.addPoint(x7,y7);
%     p7 = g.addPoint(x8,y8);
%     p8 = g.addPoint(x9,y9);
%     o = g.addPoint(0,0);
%     o2 = g.addPoint(xc,0);
% 
%     e1 = g.addSegment(p1,p2);
%     e2 = g.addSegment(p2,p3);
%     e3 = g.addSegment(p3,p4);
%     e4 = g.addArc(o,p4,p6,1);
%     e5 = g.addSegment(p6,p7);
%     e6 = g.addSegment(p7,p8);
%     e7 = g.addArc(o,p8,p1,0);   
%     e8 = g.addSegment(p2,p5);
%     e9 = g.addArc(o,p5,p7,1);
%     e10 = g.addArc(o2,p3,p5,1);
% 
%     % add loops
%     l1 = g.addLoop(e1,e8,e9,e6,e7);
%     l2 = g.addLoop(e2,e10,-e8);
%     l3 = g.addLoop(e3,e4,e5,-e9,-e10);
% 
% else
% 
%     xc = (x6^2+y6^2-x3^2)/(2*(x6-x3));   
%     p1 = g.addPoint(x1,y1);
%     p2 = g.addPoint(x2,y2);
%     p3 = g.addPoint(x3,y3);
%     p4 = g.addPoint(x4,y4);
%     p5 = g.addPoint(x5,y5);
%     p6 = g.addPoint(x6,y6);
%     p7 = g.addPoint(x7,y7);
%     p8 = g.addPoint(x8,y8);
%     p9 = g.addPoint(x9,y9);
%     o = g.addPoint(0,0);
%     o2 = g.addPoint(xc,0);
% 
%     e1 = g.addSegment(p1,p2);
%     e2 = g.addSegment(p2,p3);
%     e3 = g.addSegment(p3,p4);
%     e4 = g.addArc(o,p4,p7,1);
%     e5 = g.addSegment(p7,p8);
%     e6 = g.addSegment(p8,p9);
%     e7 = g.addArc(o,p9,p1,0);   
%     e8 = g.addSegment(p2,p5);
%     e9 = g.addArc(o,p5,p8,1);
%     e10 = g.addArc(o2,p3,p6,1);
%     e11 = g.addSegment(p6,p5);
% 
%     % add loops
%     l1 = g.addLoop(e1,e8,e9,e6,e7);
%     l2 = g.addLoop(e2,e10,e11,-e8);
%     l3 = g.addLoop(e3,e4,e5,-e9,-e11,-e10);
% 
% end
% 
% add faces
g.addFace('moving_core', l_mc);
% g.addFace(name2, l2);
% g.addFace(name3, l3);
% 
% set face colors
g.setFaceColor('moving_core',200,200,200)
% g.setFaceColor(name2,160,78,146);
% g.setFaceColor(name3,0,255,255)

% visualizations for debug
close all;
% g.setMeshMaxLength(2);
if nargin ==0, g.showSketch; end
if nargin ==0, g.showFaces; end

end
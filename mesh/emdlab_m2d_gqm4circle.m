function meshZone = emdlab_m2d_gqm4circle(x0,y0,r,meshSize)
% Generation of a structured quad mesh for a circle
% Inputs:
%   x0, y0   : circle center
%   r        : radius
%   meshSize : target element size
%
% Output:
%   meshZone : mesh object created by emdlab_m2d_qmz(elem,node)

if nargin == 0
    x0 = 0;
    y0 = 0;
    r = 5;
    meshSize = 1;
end

if meshSize <= 0
    error('meshSize must be positive.');
end
if r <= 0
    error('r must be positive.');
end

% Inner square half-size.
a = 0.45 * r;

% Divisions along each quadrant from target boundary edge length.
nArc = max(2, ceil((0.5*pi*r) / meshSize));

% Match the center block to the surrounding blocks.
n = nArc;

% Radial layers between inner square and circle.
ringThickness = r - a;
nRad = max(2, ceil(ringThickness / meshSize));

node = zeros(0,2);
elem = zeros(0,4);

tol = 1e-12;
nodeMap = containers.Map('KeyType','char','ValueType','int32');

    function id = getNode(x,y)
        xr = round(x/tol)*tol;
        yr = round(y/tol)*tol;
        key = sprintf('%.12g,%.12g', xr, yr);
        if isKey(nodeMap, key)
            id = nodeMap(key);
        else
            node(end+1,:) = [x, y];
            id = size(node,1);
            nodeMap(key) = id;
        end
    end

    function q = orientQuadCCW(q)
        p = node(q,:);
        twiceArea = sum(p(:,1) .* p([2 3 4 1],2) - p([2 3 4 1],1) .* p(:,2));
        if twiceArea < 0
            q = q([1 4 3 2]);
        end
    end

    function addGrid(X,Y)
        [nj,ni] = size(X);
        ids = zeros(nj,ni);

        for j = 1:nj
            for i = 1:ni
                ids(j,i) = getNode(X(j,i),Y(j,i));
            end
        end

        for j = 1:nj-1
            for i = 1:ni-1
                q = [ids(j,i), ids(j,i+1), ids(j+1,i+1), ids(j+1,i)];
                elem(end+1,:) = orientQuadCCW(q);
            end
        end
    end

    function p = arc(theta)
        p = [x0 + r*cos(theta), y0 + r*sin(theta)];
    end

    function addTFIBlock(Cleft,Cright,Cbot,Ctop,ns,nt)
        s = linspace(0,1,ns);
        t = linspace(0,1,nt);
        [S,T] = meshgrid(s,t);

        P00 = Cbot(0);
        P10 = Cbot(1);
        P01 = Ctop(0);
        P11 = Ctop(1);

        X = zeros(nt,ns);
        Y = zeros(nt,ns);

        for j = 1:nt
            for i = 1:ns
                ss = S(j,i);
                tt = T(j,i);

                pb = Cbot(ss);
                pt = Ctop(ss);
                pl = Cleft(tt);
                pr = Cright(tt);

                p = (1-tt)*pb + tt*pt + (1-ss)*pl + ss*pr ...
                    - ((1-ss)*(1-tt)*P00 + ss*(1-tt)*P10 + (1-ss)*tt*P01 + ss*tt*P11);

                X(j,i) = p(1);
                Y(j,i) = p(2);
            end
        end

        addGrid(X,Y);
    end

% Inner square corners
SW = [x0-a, y0-a];
SE = [x0+a, y0-a];
NE = [x0+a, y0+a];
NW = [x0-a, y0+a];

% Center block
s = linspace(0,1,n+1);
t = linspace(0,1,n+1);
[S,T] = meshgrid(s,t);

Xc = (1-S).*(1-T)*SW(1) + S.*(1-T)*SE(1) + S.*T*NE(1) + (1-S).*T*NW(1);
Yc = (1-S).*(1-T)*SW(2) + S.*(1-T)*SE(2) + S.*T*NE(2) + (1-S).*T*NW(2);

addGrid(Xc,Yc);

% East block
Cleft  = @(t) (1-t)*SE + t*NE;
Cright = @(t) arc(-pi/4 + t*(pi/2));
Cbot   = @(s) (1-s)*SE + s*arc(-pi/4);
Ctop   = @(s) (1-s)*NE + s*arc(pi/4);
addTFIBlock(Cleft,Cright,Cbot,Ctop,nRad+1,n+1);

% North block
Cleft  = @(t) (1-t)*NE + t*arc(pi/4);
Cright = @(t) (1-t)*NW + t*arc(3*pi/4);
Cbot   = @(s) (1-s)*NE + s*NW;
Ctop   = @(s) arc(pi/4 + s*(pi/2));
addTFIBlock(Cleft,Cright,Cbot,Ctop,n+1,nRad+1);

% West block
Cleft  = @(t) arc(3*pi/4 + t*(pi/2));
Cright = @(t) (1-t)*NW + t*SW;
Cbot   = @(s) (1-s)*arc(3*pi/4) + s*NW;
Ctop   = @(s) (1-s)*arc(5*pi/4) + s*SW;
addTFIBlock(Cleft,Cright,Cbot,Ctop,nRad+1,n+1);

% South block
Cleft  = @(t) (1-t)*arc(5*pi/4) + t*SW;
Cright = @(t) (1-t)*arc(7*pi/4) + t*SE;
Cbot   = @(s) arc(5*pi/4 + s*(pi/2));
Ctop   = @(s) (1-s)*SW + s*SE;
addTFIBlock(Cleft,Cright,Cbot,Ctop,n+1,nRad+1);

node = emdlab_m2d_smoothqm(node, elem, 20, 0.5);
meshZone = emdlab_m2d_qmz(elem, node);

if nargin == 0
    meshZone.showm;
end

end

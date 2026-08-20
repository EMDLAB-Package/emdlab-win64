function MeshZone = emdlab_m2d_mg0(f, v, f1, v1, maxIteration)

f1 = f;
v1 = v;

% validate the boundary
emdlab_g2d_validateFV(f,v);

% set default number of max iterations
if nargin < 5
    maxIteration = 100;
end

% length of boundary edges
el = sqrt(sum((v(f(:,1),:)-v(f(:,2),:)).^2,2));

% mid point of edges
midp = (v(f(:,1),:)+v(f(:,2),:))/2;

% edge length function
fh = scatteredInterpolant(midp(:,1),midp(:,2),el);
fh = @(p) fh(p(:,1),p(:,2));

% minimum edge length
h0 = min(el);

% Parameters of the mesh generator
dptol=.005; ttol=.05; Fscale=1.2; deltat=0.2; geps=1e-6; deps=sqrt(eps)*h0;
densityctrlfreq=30;

% 1. Create initial distribution in bounding box (equilateral triangles)
bbox = [min(v);max(v)];
[x,y] = meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
% Shift even rows
x(2:2:end,:) = x(2:2:end,:)+h0/2;
% List of node coordinates
p=[x(:),y(:)];

% 2. Remove points outside the region, apply the rejection method
% Keep only d<0 points
p = p(ext_dpoly2d(p,f1,v1)<-geps,:);
% Probability to keep point
r0 = 1./fh(p).^2;
% Rejection method
p = p(rand(size(p,1),1)<r0./max(r0),:);

% Prepend fix points
nfix = size(v,1);
p = [v; p];
% Number of points
N = size(p,1);

% 3. loop for moving nodes
% number of iterations
count = 0;
% For first iteration
pold = inf;
while count <= maxIteration
    count = count+1;    

    % 4. Retriangulation by the Delaunay algorithm
    % Any large movement?
    if max(sqrt(sum((p-pold).^2,2))/h0)>ttol

        % Save current positions
        pold = p;

        % List of triangles
        t = delaunay(p);

        % Compute centroids
        pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;

        % Keep interior triangles
        t = t(ext_dpoly2d(pmid,f1,v1)<-geps,:);

        % 5. Describe each bar by a unique pair of nodes
        % Interior bars duplicated
        bars = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];

        % Bars as node pairs
        bars = unique(sort(bars,2),'rows');

    end

    % 6. Move mesh points based on bar lengths L and forces F
    % List of bar vectors
    barvec = p(bars(:,1),:)-p(bars(:,2),:);
    % L = Bar lengths
    L = sqrt(sum(barvec.^2,2));
    hbars=fh((p(bars(:,1),:)+p(bars(:,2),:))/2);
    % L0 = Desired lengths
    L0 = hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));
    % Density control - remove points that are too close
    if mod(count,densityctrlfreq)==0 && any(L0>2*L)
        p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nfix),:)=[];
        N=size(p,1); pold=inf;
        continue;
    end
    % Bar forces (scalars)
    F = max(L0-L,0);
    % Bar forces (x,y components)
    Fvec=F./L*[1,1].*barvec;
    Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    % Force = 0 at fixed points
    Ftot(1:size(v,1),:)=0;
    % Update node positions
    p=p+deltat*Ftot;
    
    % 7. Bring outside points back to the boundary
    % Find points outside (d>0)
    d=ext_dpoly2d(p,f1,v1); ix=d>geps;
    % Numerical gradient
    dgradx=(ext_dpoly2d([p(ix,1)+deps,p(ix,2)],f1,v1)-d(ix))/deps;
    dgrady=(ext_dpoly2d([p(ix,1),p(ix,2)+deps],f1,v1)-d(ix))/deps;
    dgrad2=dgradx.^2+dgrady.^2;
    % Project
    p(ix,:)=p(ix,:)-[d(ix).*dgradx./dgrad2,d(ix).*dgrady./dgrad2];

    % 8. Termination criterion: All interior nodes move less than dptol (scaled)
    if max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0)<dptol
        break
    end

    % minimal mesh
    if size(p,1) == size(v,1)
        break;
    end    

end

% 9. Fix mesh
p = p(ext_dpoly2d(p,f1,v1)<-geps,:);
p = [v;p];
t = delaunayTriangulation(p,f);
t = t.ConnectivityList(ext_dpoly2d(t.incenter,f1,v1)<-geps,:);

MeshZone = emdlab_m2d_tmz(t,p);
MeshZone.moveNodes;

end

function [quads, triLeft] = emdlab_m2d_pair_triangles(tcl, pts)
% Pair adjacent triangles into quads using a geometric quality score.
%
% Inputs
%   tcl : NT x 3 triangle connectivity
%   pts : NP x 2 or NP x 3 point coordinates
%
% Outputs
%   quads   : NQ x 4 quad connectivity
%   triLeft : remaining unpaired triangles

if size(pts,2) > 2
    pts = pts(:,1:2);
end

nt = size(tcl,1);

% Build undirected triangle edges
E = [tcl(:,[1 2]);
     tcl(:,[2 3]);
     tcl(:,[3 1])];

triId = [(1:nt)'; (1:nt)'; (1:nt)'];
Es = sort(E,2);

[Esu,~,ic] = unique(Es,'rows');
counts = accumarray(ic,1);

% Only interior edges can join two triangles
intEdgeIds = find(counts == 2);

% Candidate rows: [triA triB score q1 q2 q3 q4]
candidates = zeros(0,7);

for ue = intEdgeIds'
    idx = find(ic == ue);
    tA = triId(idx(1));
    tB = triId(idx(2));

    e = Esu(ue,:);
    i = e(1);
    j = e(2);

    aNodes = tcl(tA,:);
    bNodes = tcl(tB,:);

    k = aNodes(aNodes ~= i & aNodes ~= j);
    l = bNodes(bNodes ~= i & bNodes ~= j);

    if isempty(k) || isempty(l) || k == l
        continue
    end

    % Two possible cyclic orderings
    q1 = [k i l j];
    q2 = [k j l i];

    [ok1,s1] = quadScore(q1, pts);
    [ok2,s2] = quadScore(q2, pts);

    if ok1 || ok2
        if ok1 && (~ok2 || s1 >= s2)
            candidates(end+1,:) = [tA tB s1 q1]; %#ok<AGROW>
        else
            candidates(end+1,:) = [tA tB s2 q2]; %#ok<AGROW>
        end
    end
end

if isempty(candidates)
    quads = zeros(0,4);
    triLeft = tcl;
    return
end

% Greedy maximum-weight selection
[~,ord] = sort(candidates(:,3),'descend');
cand = candidates(ord,:);

used = false(nt,1);
quads = zeros(0,4);

for r = 1:size(cand,1)
    tA = cand(r,1);
    tB = cand(r,2);

    if ~used(tA) && ~used(tB)
        used(tA) = true;
        used(tB) = true;
        quads(end+1,:) = cand(r,4:7); %#ok<AGROW>
    end
end

triLeft = tcl(~used,:);

end

function [ok,score] = quadScore(q, pts)
p = pts(q,:);

% Reject repeated nodes
if numel(unique(q)) < 4
    ok = false;
    score = -inf;
    return
end

% Reject self-intersecting or degenerate quads
if ~isSimpleQuad(p)
    ok = false;
    score = -inf;
    return
end

% Enforce convexity for safer FE / meshing behavior
if ~isConvexQuad(p)
    ok = false;
    score = -inf;
    return
end

A = polyarea(p(:,1), p(:,2));
if A <= 0
    ok = false;
    score = -inf;
    return
end

theta = quadAngles(p) * 180/pi;
L = sqrt(sum((p([2 3 4 1],:) - p).^2, 2));

% Hard quality filters
if any(theta < 20) || any(theta > 160)
    ok = false;
    score = -inf;
    return
end

lr = max(L) / max(min(L), eps);
if lr > 4.0
    ok = false;
    score = -inf;
    return
end

d1 = norm(p(3,:) - p(1,:));
d2 = norm(p(4,:) - p(2,:));
dr = max(d1,d2) / max(min(d1,d2), eps);

% Lower penalties are better
anglePenalty = sum((theta - 90).^2);
lenPenalty   = (lr - 1)^2;
diagPenalty  = (dr - 1)^2;

% Higher score is better
score = 1 / (1e-12 + anglePenalty + 2*lenPenalty + 0.5*diagPenalty);
ok = true;
end

function tf = isSimpleQuad(p)
% Check for nonzero edge lengths
e = p([2 3 4 1],:) - p;
if any(sqrt(sum(e.^2,2)) <= 1e-14)
    tf = false;
    return
end

% A quad is simple iff its non-adjacent edges do not intersect
tf = true;

if segmentsIntersect(p(1,:), p(2,:), p(3,:), p(4,:))
    tf = false;
    return
end

if segmentsIntersect(p(2,:), p(3,:), p(4,:), p(1,:))
    tf = false;
    return
end
end

function tf = isConvexQuad(p)
% Signed z-components of successive cross products must have same sign
c = zeros(4,1);
for k = 1:4
    a = p(mod(k,4)+1,:) - p(k,:);
    b = p(mod(k+1,4)+1,:) - p(mod(k,4)+1,:);
    c(k) = cross2d(a,b);
end

tol = 1e-12;
if any(abs(c) <= tol)
    tf = false;
    return
end

tf = all(c > 0) || all(c < 0);
end

function theta = quadAngles(p)
theta = zeros(4,1);
for k = 1:4
    km = mod(k-2,4) + 1;
    kp = mod(k,4) + 1;

    v1 = p(km,:) - p(k,:);
    v2 = p(kp,:) - p(k,:);

    n1 = norm(v1);
    n2 = norm(v2);
    if n1 <= eps || n2 <= eps
        theta(k) = pi;
        continue
    end

    c = dot(v1,v2) / (n1*n2);
    c = max(-1,min(1,c));
    theta(k) = acos(c);
end
end

function tf = segmentsIntersect(a,b,c,d)
% Proper intersection test for two closed segments.
o1 = orient2d(a,b,c);
o2 = orient2d(a,b,d);
o3 = orient2d(c,d,a);
o4 = orient2d(c,d,b);

tol = 1e-12;

% Proper crossing
if ((o1 > tol && o2 < -tol) || (o1 < -tol && o2 > tol)) && ...
   ((o3 > tol && o4 < -tol) || (o3 < -tol && o4 > tol))
    tf = true;
    return
end

% Collinear / touching cases
tf = (abs(o1) <= tol && onSegment(a,b,c)) || ...
     (abs(o2) <= tol && onSegment(a,b,d)) || ...
     (abs(o3) <= tol && onSegment(c,d,a)) || ...
     (abs(o4) <= tol && onSegment(c,d,b));
end

function tf = onSegment(a,b,p)
tol = 1e-12;
tf = p(1) >= min(a(1),b(1)) - tol && p(1) <= max(a(1),b(1)) + tol && ...
     p(2) >= min(a(2),b(2)) - tol && p(2) <= max(a(2),b(2)) + tol;
end

function z = orient2d(a,b,c)
z = cross2d(b-a, c-a);
end

function z = cross2d(u,v)
z = u(1)*v(2) - u(2)*v(1);
end

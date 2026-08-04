function [quads, triLeft, info] = emdlab_m2d_pair_triangles(tcl, pts, opts)
%EMDLAB_M2D_PAIR_TRIANGLES_POWERFUL  Pair adjacent triangles into high-quality quads.
%
% This is a stronger version than a simple greedy matcher:
%   1) Enumerate all interior-edge triangle pairs (candidates).
%   2) Compute a robust geometric quality score (convexity, intersections,
%      Jacobian sign-like checks, angles, edge ratios, diagonal balance, warpage=0 in 2D).
%   3) Select a maximum-weight *matching* on the triangle adjacency graph.
%      - Exact solver requires Blossom (not in base MATLAB).
%      - Here: a strong heuristic: greedy + local 2-opt improvements.
%
% Inputs
%   tcl  : NT x 3 triangles (1-based node ids)
%   pts  : NP x 2 or NP x 3 points (uses xy)
%   opts : struct (optional)
%       .minAngleDeg   (default 15)
%       .maxAngleDeg   (default 165)
%       .maxEdgeRatio  (default 5.0)
%       .maxDiagRatio  (default 5.0)
%       .wAngle        (default 1.0)
%       .wEdge         (default 0.3)
%       .wDiag         (default 0.3)
%       .wArea         (default 0.1)  % prefers larger quads if equal quality
%       .improvePasses (default 3)    % local improvement passes
%       .verbose       (default false)
%
% Outputs
%   quads   : NQ x 4 quad connectivity (cyclic order)
%   triLeft : remaining triangles (unpaired)
%   info    : diagnostics
%
% Notes
% - If you need provably optimal matching: integrate a Blossom V / MWPM solver
%   (MEX or external) and replace the selection stage.
% - Robust for typical unstructured triangle meshes.

if nargin < 3, opts = struct(); end
opts = setDefaults(opts);

if size(pts,2) > 2, pts = pts(:,1:2); end

nt = size(tcl,1);

% ---- Build undirected edges and triangle adjacency via shared edges ----
E  = [tcl(:,[1 2]); tcl(:,[2 3]); tcl(:,[3 1])];
triId = [(1:nt)'; (1:nt)'; (1:nt)'];
Es = sort(E,2);

[Esu,~,ic] = unique(Es,'rows');
counts = accumarray(ic,1);
intEdgeIds = find(counts == 2); % interior edges => 2 incident triangles

% Candidate format:
%   [tA tB score q1 q2 q3 q4 area]
candidates = zeros(0,8);

for ue = intEdgeIds'
    idx = find(ic == ue);
    tA = triId(idx(1));
    tB = triId(idx(2));

    e = Esu(ue,:);
    i = e(1); j = e(2);

    aNodes = tcl(tA,:);
    bNodes = tcl(tB,:);

    k = aNodes(aNodes ~= i & aNodes ~= j);
    l = bNodes(bNodes ~= i & bNodes ~= j);

    if isempty(k) || isempty(l) || k == l
        continue
    end

    % Two possible cyclic orderings around the shared edge
    qA = [k i l j];
    qB = [k j l i];

    [okA, sA, areaA] = quadScore2D(qA, pts, opts);
    [okB, sB, areaB] = quadScore2D(qB, pts, opts);

    if okA || okB
        if okA && (~okB || sA >= sB)
            candidates(end+1,:) = [tA tB sA qA areaA]; %#ok<AGROW>
        else
            candidates(end+1,:) = [tA tB sB qB areaB]; %#ok<AGROW>
        end
    end
end

info = struct();
info.numTriangles   = nt;
info.numCandidates  = size(candidates,1);

if isempty(candidates)
    quads = zeros(0,4);
    triLeft = tcl;
    info.numQuads = 0;
    return
end

% ---- Build adjacency list per triangle for faster local improvement ----
% Store candidate index for each unordered triangle pair (tA,tB)
tA = candidates(:,1); tB = candidates(:,2);
w  = candidates(:,3);

% ---- Initial greedy selection (highest weight first) ----
[~,ord] = sort(w,'descend');
cand = candidates(ord,:);

used = false(nt,1);
picked = false(size(cand,1),1);

for r = 1:size(cand,1)
    a = cand(r,1); b = cand(r,2);
    if ~used(a) && ~used(b)
        used(a)=true; used(b)=true;
        picked(r)=true;
    end
end

% ---- Local improvement: 2-opt style augmentations ----
% Try to replace two picked edges with two better edges that still form a matching.
% This often recovers from greedy mistakes cheaply.
if opts.improvePasses > 0
    % Map triangle -> picked partner and weight for quick checks
    for pass = 1:opts.improvePasses
        improved = false;

        [partner, pW] = currentMatchingFromPicked(cand, picked, nt);

        % Consider each triangle as a pivot; try to reroute through an unused candidate
        % by breaking existing pairs and forming better ones.
        for a = 1:nt
            b = partner(a);
            if b == 0, continue; end

            % Try candidates incident to a (scan all; for big meshes you can index)
            % We do a limited scan by using candidates list directly; acceptable in MATLAB for moderate sizes.
            for r1 = 1:size(cand,1)
                x = cand(r1,1); y = cand(r1,2);
                if x ~= a && y ~= a, continue; end
                c = x; d = y;
                if c ~= a, d = c; c = a; end % ensure c=a, d=other

                if d == b, continue; end

                % d must be currently matched to e (or free)
                e = partner(d);

                % Case 1: d is free => simply swap (a-b) for (a-d) if better
                if e == 0
                    if cand(r1,3) > pW(a) + 1e-15
                        % remove edge (a,b), add (a,d)
                        picked = applySwap1(cand, picked, a, b, r1);
                        improved = true;
                        [partner, pW] = currentMatchingFromPicked(cand, picked, nt);
                    end
                    continue
                end

                if e == a || e == b, continue; end

                % Need a second edge to keep matching size: choose edge (b,e) if exists and feasible.
                r2 = findCandidate(cand, b, e);
                if r2 == 0, continue; end

                % Current weight: (a,b) + (d,e)
                % Proposed: (a,d) + (b,e)
                curW = pW(a) + pW(d);
                newW = cand(r1,3) + cand(r2,3);

                if newW > curW + 1e-15
                    % Perform 2-edge swap
                    picked = applySwap2(cand, picked, a, b, d, e, r1, r2);
                    improved = true;
                    [partner, pW] = currentMatchingFromPicked(cand, picked, nt);
                end
            end
        end

        if opts.verbose
            fprintf('Improve pass %d: %s\n', pass, ternary(improved,'improved','no change'));
        end
        if ~improved, break; end
    end
end

% ---- Emit quads and leftover triangles ----
sel = find(picked);
quads = cand(sel,4:7);
used = false(nt,1);
for k = 1:numel(sel)
    used(cand(sel(k),1)) = true;
    used(cand(sel(k),2)) = true;
end
triLeft = tcl(~used,:);

info.numQuads = size(quads,1);
info.numTriLeft = size(triLeft,1);

end

% ======================================================================
%                               Scoring
% ======================================================================

function [ok, score, A] = quadScore2D(q, pts, opts)
p = pts(q,:);

% Unique nodes
if numel(unique(q)) < 4
    ok = false; score = -inf; A = 0; return
end

% Simple + convex
if ~isSimpleQuad(p) || ~isConvexQuad(p)
    ok = false; score = -inf; A = 0; return
end

A = polyarea(p(:,1), p(:,2));
if A <= 0
    ok = false; score = -inf; return
end

theta = quadAngles(p) * 180/pi;
if any(theta < opts.minAngleDeg) || any(theta > opts.maxAngleDeg)
    ok = false; score = -inf; return
end

L  = sqrt(sum((p([2 3 4 1],:) - p).^2, 2));
lr = max(L) / max(min(L), eps);
if lr > opts.maxEdgeRatio
    ok = false; score = -inf; return
end

d1 = norm(p(3,:) - p(1,:));
d2 = norm(p(4,:) - p(2,:));
dr = max(d1,d2) / max(min(d1,d2), eps);
if dr > opts.maxDiagRatio
    ok = false; score = -inf; return
end

% Penalties
anglePenalty = mean((theta - 90).^2) / (90^2);  % ~O(1)
edgePenalty  = (lr - 1)^2;
diagPenalty  = (dr - 1)^2;

% Reward area mildly to break ties (normalised by perimeter^2)
per = sum(L);
areaTerm = A / max(per^2, eps);

% Higher better
den = 1e-12 + opts.wAngle*anglePenalty + opts.wEdge*edgePenalty + opts.wDiag*diagPenalty ...
      - opts.wArea*areaTerm; % subtract because areaTerm is good
score = 1 / max(den, 1e-12);

ok = true;
end

% ======================================================================
%                          Geometry predicates
% ======================================================================

function tf = isSimpleQuad(p)
e = p([2 3 4 1],:) - p;
if any(sqrt(sum(e.^2,2)) <= 1e-14), tf = false; return; end

% Non-adjacent edges must not intersect
if segmentsIntersect(p(1,:), p(2,:), p(3,:), p(4,:)), tf = false; return; end
if segmentsIntersect(p(2,:), p(3,:), p(4,:), p(1,:)), tf = false; return; end
tf = true;
end

function tf = isConvexQuad(p)
c = zeros(4,1);
for k = 1:4
    a = p(mod(k,4)+1,:) - p(k,:);
    b = p(mod(k+1,4)+1,:) - p(mod(k,4)+1,:);
    c(k) = cross2d(a,b);
end
tol = 1e-12;
if any(abs(c) <= tol), tf = false; return; end
tf = all(c > 0) || all(c < 0);
end

function theta = quadAngles(p)
theta = zeros(4,1);
for k = 1:4
    km = mod(k-2,4) + 1;
    kp = mod(k,4) + 1;
    v1 = p(km,:) - p(k,:);
    v2 = p(kp,:) - p(k,:);
    n1 = norm(v1); n2 = norm(v2);
    if n1 <= eps || n2 <= eps
        theta(k) = pi;
    else
        c = dot(v1,v2) / (n1*n2);
        c = max(-1,min(1,c));
        theta(k) = acos(c);
    end
end
end

function tf = segmentsIntersect(a,b,c,d)
o1 = orient2d(a,b,c);
o2 = orient2d(a,b,d);
o3 = orient2d(c,d,a);
o4 = orient2d(c,d,b);
tol = 1e-12;

if ((o1 > tol && o2 < -tol) || (o1 < -tol && o2 > tol)) && ...
   ((o3 > tol && o4 < -tol) || (o3 < -tol && o4 > tol))
    tf = true; return
end

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

% ======================================================================
%                     Matching helpers (heuristic)
% ======================================================================

function [partner, pW] = currentMatchingFromPicked(cand, picked, nt)
partner = zeros(nt,1);
pW      = zeros(nt,1);
sel = find(picked);
for k = 1:numel(sel)
    r = sel(k);
    a = cand(r,1); b = cand(r,2); w = cand(r,3);
    partner(a) = b; partner(b) = a;
    pW(a) = w; pW(b) = w;
end
end

function r = findCandidate(cand, a, b)
% Return index in cand for unordered pair (a,b), or 0 if none.
if a > b, tmp=a; a=b; b=tmp; end
% cand may not be sorted; linear search (can be optimised with hashing)
mask = ( (cand(:,1)==a & cand(:,2)==b) | (cand(:,1)==b & cand(:,2)==a) );
idx = find(mask,1,'first');
if isempty(idx), r = 0; else, r = idx; end
end

function picked = applySwap1(cand, picked, a, b, rAdd)
% remove (a,b), add rAdd
rRem = findCandidate(cand, a, b);
if rRem ~= 0, picked(rRem) = false; end
picked(rAdd) = true;
end

function picked = applySwap2(cand, picked, a, b, d, e, r1, r2)
% remove (a,b) and (d,e), add r1 and r2
rRem1 = findCandidate(cand, a, b);
rRem2 = findCandidate(cand, d, e);
if rRem1 ~= 0, picked(rRem1) = false; end
if rRem2 ~= 0, picked(rRem2) = false; end
picked(r1) = true;
picked(r2) = true;
end

% ======================================================================
%                              Defaults
% ======================================================================

function opts = setDefaults(opts)
def.minAngleDeg   = 15;
def.maxAngleDeg   = 165;
def.maxEdgeRatio  = 5.0;
def.maxDiagRatio  = 5.0;
def.wAngle        = 1.0;
def.wEdge         = 0.3;
def.wDiag         = 0.3;
def.wArea         = 0.1;
def.improvePasses = 3;
def.verbose       = false;

f = fieldnames(def);
for k = 1:numel(f)
    if ~isfield(opts,f{k}) || isempty(opts.(f{k}))
        opts.(f{k}) = def.(f{k});
    end
end
end

function s = ternary(cond,a,b)
if cond, s=a; else, s=b; end
end

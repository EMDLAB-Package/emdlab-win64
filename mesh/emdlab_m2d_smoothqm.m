function node = emdlab_m2d_smoothqm(node, elem, nIter, alpha)
if nargin < 3
    nIter = 20;
end
if nargin < 4
    alpha = 0.4;
end

nNode = size(node,1);

% Build node adjacency
nbr = cell(nNode,1);
for e = 1:size(elem,1)
    q = elem(e,:);
    edges = [q([1 2]); q([2 3]); q([3 4]); q([4 1])];
    for k = 1:4
        i = edges(k,1);
        j = edges(k,2);
        nbr{i}(end+1) = j;
        nbr{j}(end+1) = i;
    end
end

for i = 1:nNode
    nbr{i} = unique(nbr{i});
end

% Detect boundary nodes from edges used by only one element
edgeList = zeros(size(elem,1)*4, 2);
row = 0;
for e = 1:size(elem,1)
    q = elem(e,:);
    localEdges = [q([1 2]); q([2 3]); q([3 4]); q([4 1])];
    for k = 1:4
        row = row + 1;
        edgeList(row,:) = sort(localEdges(k,:));
    end
end

[uniqueEdges, ~, ic] = unique(edgeList, 'rows');
counts = accumarray(ic, 1);
bndEdges = uniqueEdges(counts == 1, :);
isBoundary = false(nNode,1);
isBoundary(unique(bndEdges(:))) = true;

% Laplacian smoothing with fixed boundary
for it = 1:nIter
    newNode = node;

    for i = 1:nNode
        if isBoundary(i)
            continue;
        end

        nei = nbr{i};
        if isempty(nei)
            continue;
        end

        target = mean(node(nei,:), 1);
        candidate = node(i,:) + alpha * (target - node(i,:));

        % Optional: accept only if local elements stay valid
        if isValidMove(i, candidate, node, elem)
            newNode(i,:) = candidate;
        end
    end

    node = newNode;
end
end

function ok = isValidMove(iNode, candidate, node, elem)
ok = true;
attached = any(elem == iNode, 2);
elemIds = find(attached);

testNode = node;
testNode(iNode,:) = candidate;

for k = 1:numel(elemIds)
    q = elem(elemIds(k),:);
    p = testNode(q,:);
    area2 = sum(p(:,1).*p([2 3 4 1],2) - p([2 3 4 1],1).*p(:,2));
    if area2 <= 1e-12
        ok = false;
        return;
    end
end
end

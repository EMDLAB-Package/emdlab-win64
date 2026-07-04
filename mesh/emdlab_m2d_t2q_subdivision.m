function [cells_quad,points_quad] = emdlab_m2d_t2q_subdivision(triangles,points)
    % Inputs:
    %   points    : Nn x 2 matrix of vertex coordinates [x, y]
    %   triangles : Ne x 3 matrix of triangle connectivity (indices)
    %
    % Outputs:
    %   points_quad : New vertex coordinates (Nn_new x 2)
    %   cells_quad  : Quadrilateral connectivity matrix (3*Ne x 4)

    Nn = size(points, 1);
    Ne = size(triangles, 1);

    % 1. Extract unique edges and identify midpoints without duplication
    edges = [triangles(:,[1 2]); triangles(:,[2 3]); triangles(:,[3 1])];
    edges_sorted = sort(edges, 2);
    [unique_edges, ~, edge_map] = unique(edges_sorted, 'rows');
    num_edges = size(unique_edges, 1);

    % Compute midpoint coordinates
    midpoints = (points(unique_edges(:,1), :) + points(unique_edges(:,2), :)) / 2;

    % 2. Compute triangle centroids
    centroids = zeros(Ne, 2);
    for i = 1:Ne
        centroids(i, :) = mean(points(triangles(i,:), :), 1);
    end

    % 3. Assemble new points matrix
    % Order: [Original Points; Midpoint Points; Centroid Points]
    points_quad = [points; midpoints; centroids];

    % Index offsets to reference midpoints and centroids
    idx_midpoint_offset = Nn;
    idx_centroid_offset = Nn + num_edges;

    % 4. Build the new quadrilateral connectivity matrix
    cells_quad = zeros(Ne * 3, 4);

    % Map back the original edge mappings to retrieve local edge midpoints
    map_e12 = edge_map(1:Ne);
    map_e23 = edge_map(Ne+1 : 2*Ne);
    map_e31 = edge_map(2*Ne+1 : end);

    for i = 1:Ne
        v1 = triangles(i, 1);
        v2 = triangles(i, 2);
        v3 = triangles(i, 3);

        % Get global indices of the three edge midpoints
        m12 = map_e12(i) + idx_midpoint_offset;
        m23 = map_e23(i) + idx_midpoint_offset;
        m31 = map_e31(i) + idx_midpoint_offset;

        % Get global index of this triangle's centroid
        c = i + idx_centroid_offset;

        % Subdivide triangle into 3 quadrilaterals
        cells_quad((i-1)*3 + 1, :) = [v1, m12, c, m31];
        cells_quad((i-1)*3 + 2, :) = [m12, v2, m23, c];
        cells_quad((i-1)*3 + 3, :) = [c, m23, v3, m31];
    end
end

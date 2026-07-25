% EMDLAB: Electrical Machines Design Laboratory
% Prism mesh zone (3D element)

classdef emdlab_m3d_pmz < handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_m3d_xmz
    
    properties (SetAccess = protected)

        % tfacets
        tfacets(3, :) double;

        % qfacets
        qfacets(4, :) double;

    end
    
    properties (Dependent = true)

        % number of triangular facets
        Ntf(1, 1) double;

        % number of quadrilateral facets
        Nqf(1, 1) double;

    end
    
    methods
        %% Constructor and Destructor
        function obj = emdlab_m3d_pmz(cl, nodes)

            if nargin < 2, error('Not enough input arguments.'); end
            if nargin > 2, error('Too many input arguments.'); end

            obj.cl = cl;
            obj.nodes = nodes;

        end
        
        function y = get.Ntf(obj)
            y = size(obj.qfacets, 2);
        end
        
        function y = get.Nqf(obj)
            y = size(obj.tfacets, 2);
        end

        function setData(obj)

            % check if already data is set
            if obj.isDataSet, return; end

            % tetrahedral facets
            f1 = obj.cl(:,[1,2,3]);
            f2 = obj.cl(:,[2,4,3]);
            f3 = obj.cl(:,[3,4,1]);
            f4 = obj.cl(:,[1,4,2]);

            % sorting for lower index
            [f1,s1] = sort(f1,2);
            [f2,s2] = sort(f2,2);
            [f3,s3] = sort(f3,2);
            [f4,s4] = sort(f4,2);

            % specefying changed facet index
            s1 = ((s1(:,1)==1)&(s1(:,2)==3))|...
                ((s1(:,1)==3)&(s1(:,2)==2))|...
                ((s1(:,1)==2)&(s1(:,2)==1));
            s2 = ((s2(:,1)==1)&(s2(:,2)==3))|...
                ((s2(:,1)==3)&(s2(:,2)==2))|...
                ((s2(:,1)==2)&(s2(:,2)==1));
            s3 = ((s3(:,1)==1)&(s3(:,2)==3))|...
                ((s3(:,1)==3)&(s3(:,2)==2))|...
                ((s3(:,1)==2)&(s3(:,2)==1));
            s4 = ((s4(:,1)==1)&(s4(:,2)==3))|...
                ((s4(:,1)==3)&(s4(:,2)==2))|...
                ((s4(:,1)==2)&(s4(:,2)==1));
            
            % unification of facets
            [obj.facets,~,ic] = unique([f1;f2;f3;f4],'rows');

            % getting number of elements
            ne = obj.Ne;

            % getting index of facets corresponding to each elements
            f1 = ic(1:ne);
            f2 = ic(1+ne:2*ne);
            f3 = ic(1+2*ne:3*ne);
            f4 = ic(1+3*ne:4*ne);

            % specefying boundary facets
            obj.bfacets = sparse([f1,f2,f3,f4],ones(4*ne,1),ones(4*ne,1));
            obj.bfacets = full(obj.bfacets == 1);
            
            % specefying trace direction
            f1(s1) = -f1(s1);
            f2(s2) = -f2(s2);
            f3(s3) = -f3(s3);
            f4(s4) = -f4(s4);
            
            % element matrix
            obj.elements = [f1,f2,f3,f4];
            
            % evaluation of area of each elements
            obj.calculateVolumeOfElements;
            obj.calculateMeshZoneVolume;
            
            % change states
            obj.isDataSet = true;

        end
        
        function setData_old(obj)
                        
            % check if already data is set
            if obj.isDataSet, return; end

            % tetrahedral facets
            e1 = obj.cl([1, 2], :);
            e2 = obj.cl([2, 3], :);
            e3 = obj.cl([3, 1], :);
            e4 = obj.cl([4, 5], :);
            e5 = obj.cl([5, 6], :);
            e6 = obj.cl([6, 4], :);
            e7 = obj.cl([1, 4], :);
            e8 = obj.cl([2, 5], :);
            e9 = obj.cl([3, 6], :);            
            
            % sorting for lower index
            [e1, s1] = sort(e1);
            [e2, s2] = sort(e2);
            [e3, s3] = sort(e3);
            [e4, s4] = sort(e4);
            [e5, s5] = sort(e5);
            [e6, s6] = sort(e6);
            [e7, s7] = sort(e7);
            [e8, s8] = sort(e8);
            [e9, s9] = sort(e9);
            
            % specefying changed edge index
            s1 = s1(1,:) == 2;
            s2 = s2(1,:) == 2;
            s3 = s3(1,:) == 2;
            s4 = s4(1,:) == 2;
            s5 = s5(1,:) == 2;
            s6 = s6(1,:) == 2;
            s7 = s7(1,:) == 2;
            s8 = s8(1,:) == 2;
            s9 = s9(1,:) == 2;
            
            % unification of edges
            [tmp, ~, ic] = unique([e1, e2, e3, e4, e5, e6, e7, e8, e9]', 'rows');
            
            obj.edges = tmp';
            % getting number of elements
            ne = obj.Ne;
            % getting index of edge corresponding to each elements
            e1 = ic(1:ne);
            e2 = ic(1 + ne:2 * ne);
            e3 = ic(1 + 2 * ne:3 * ne);
            e4 = ic(1 + 3 * ne:4 * ne);
            e5 = ic(1 + 4 * ne:5 * ne);
            e6 = ic(1 + 5 * ne:6 * ne);
            e7 = ic(1 + 6 * ne:7 * ne);
            e8 = ic(1 + 7 * ne:8 * ne);
            e9 = ic(1 + 8 * ne:9 * ne);
            
            e1(s1) = -e1(s1);
            e2(s2) = -e2(s2);
            e3(s3) = -e3(s3);
            e4(s4) = -e4(s4);
            e5(s5) = -e5(s5);
            e6(s6) = -e6(s6);
            e7(s7) = -e7(s7);
            e8(s8) = -e8(s8);
            e9(s9) = -e9(s9);
            
            obj.tfacets = [e1,e2,e3;e4,e5,e6]';
            obj.qfacets = [e1,e8,e4,e7;e2,e9,e5,e8;e3,e7,e6,e9]';
            
        end

    end
    
end

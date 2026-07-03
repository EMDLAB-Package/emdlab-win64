function [facets, nodes] = emdlab_g3d_readBinarySTL(FileDir)

f = fopen(FileDir,'r');

fread(f,80,'*char');

Nt = fread(f,1,'uint32');

n = zeros(3*Nt,1);
p = zeros(9*Nt,1);

for i = 1:Nt
    n(3*i-2:3*i) = fread(f,3,'single');
    p(9*i-8:9*i-6) = fread(f,3,'single');
    p(9*i-5:9*i-3) = fread(f,3,'single');
    p(9*i-2:9*i) = fread(f,3,'single');
    fread(f,1,'uint16');
end

fclose(f);

p = reshape(p,3,[]);
p = p';

facets = 1:3*Nt;
facets = reshape(facets,3,[]);
facets = facets';

[nodes,~,index] = uniquetol(p,1e-4,'ByRows',true);

facets = index(facets);

end
function mz = emdlab_m2d_gqm4rectangle(x0,y0,W,H,Nx,Ny)

if nargin == 0
    x0 = 0;
    y0 = 0;
    W = 1;
    H = 1;
    Nx = 10;
    Ny = 10;
end

% mesh generation
x = x0 + linspace(0,W,Nx+1);
y = y0 + linspace(0,H,Ny+1);

[x,y] = ndgrid(x,y);
x = x(:);
y = y(:);
pts = [x,y];

% connectivity list
index = 0;
cl = zeros(Nx*Ny,4);
for j = 1:Ny
    for i = 1:Nx
        index = index + 1;
        indexE = (j-1)*(Nx+1) + i;
        cl(index,:) = [indexE,indexE+1,indexE+Nx+2,indexE+Nx+1];
    end
end

mz = emdlab_m2d_qmz(cl,pts);

if nargin == 0, mz.showm; end

end
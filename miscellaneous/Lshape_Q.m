
p = [0,0
    2,0
    2,2
    1,2
    1,1
    0,1];

q = [1 2 5 6
    2 3 4 5];

plotQmesh(p,q)

%%
[p,q] = strefineQ(p,q);
[p,q] = strefineQ(p,q);
[p,q] = strefineQ(p,q);
plotQmesh(p,q)
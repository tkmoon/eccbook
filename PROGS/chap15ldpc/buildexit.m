function z=buildexit(ctob,btoc)
% function buildexit(ctob,btoc)
%
% given the mutual information ctob and btoc as a function of iteration,
% stack these up according to the iteration number to 
% form the exit chart


xpath=[0 0];
ypath=[0 btoc(1)];
for(i=1:(length(ctob)-1))
    xpath=[xpath ctob(i) ctob(i)];
    ypath=[ypath btoc(i) btoc(i+1)];
end
xpath=[xpath ctob(end)];
ypath=[ypath btoc(end)];


z=[ xpath ; ypath ];
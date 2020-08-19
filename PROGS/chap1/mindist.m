function dmin = mindist(G)
% function dmin = mindist(G)
%
% Exhuastively compute the minimum distance of a linear code described
% by the generator matrix G.
% G is an NxK matrix  (watch the convention!)
%
% Note: since this is an exhaustive test, it is probably not practical
% for codes with K>15
%


[N,K] = size(G);
dmin = N;
m = zeros(K,1);
for i=1:(2^K-1)
  m = c2b(i,K)';
  c = mod(G*m,2);
c
  d = sum(c);
  if(d < dmin)
	dmin = d;
	mmin = m;
	cmin = c;
  end
end

dmin
mmin
cmin

  
  
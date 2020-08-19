function Km = computekm(n,k,m)
% function Km = computekm(n,k,m)
% 
% Compute Km, the number of symbols correct, for a 
% (n,k) Reed-Solomon with the Guruswami(m) decoding algorithm
%

v = k-1;
if(m==0)
  Km = ceil((n+v+1)/2);
  return
end
C = n*nchoosek(m+1,2);
k1 = 0;
while(1)
  k1 = k1+1;
  l = (m*k1)-1;
  lk = floor(l/v);
  n1 = (l+1)*(lk+1) - v/2*lk*(lk+1);
  if(n1 > C) break; end
end
Km = k1;  
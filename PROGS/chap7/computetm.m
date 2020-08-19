function tm = computetm(n,k,m)
% function tm = computetm(n,k,m)
% 
% Compute tm, the error correction capability, for an
% (n,k) Reed-Solomon with the Guruswami(m) decoding algorithm
%

v = k-1;
if(m==0)
  Km = ceil((n+v+1)/2);
  tm = n-Km;
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
tm = n-Km;
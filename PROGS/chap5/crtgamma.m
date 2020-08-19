function gamma = crtgamma(m)
% function gamma = crtgamma(m)
%
% Compute the gammas for the CRT

r = length(m);
mp = prod(m);
for i=1:r
  [g,b,y1] = gcdint2(mp/m(i),m(i));
  gamma(i) = mp/m(i)*b;
end
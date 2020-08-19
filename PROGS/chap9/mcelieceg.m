function g = mcelieceg(u)
%
% compute the g function used by Mceliece bound

g = h2(0.5*(1-sqrt(1-u)));

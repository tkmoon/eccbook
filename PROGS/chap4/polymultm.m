function p = polymult(a,b,m)
% multipoly the polynomials p=a*b modulo m
% The polynomials are represented with the HIGHEST coefficient first
% e.g., x^2 + 3x + 4 -- [1 3 4],   x^2 + 3x --- [ 1 3 0]

if(all(a==0) | all(b==0))
  p = 0;
else
  k = length(b) - length(a);
  p = mod(conv(a,b),m);
end
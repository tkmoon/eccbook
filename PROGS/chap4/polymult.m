function p = polymult(a,b)
% multipoly the polynomials p=a*b
% The polynomials are represented with the HIGHEST coefficient first
% e.g., x^2 + 3x + 4 -- [1 3 4],   x^2 + 3x --- [ 1 3 0]

if(all(a==0) | all(b==0))
  p = 0;
else
  k = length(b) - length(a);
  p = conv(a,b);
end
function p = polyadd(a,b)
% add the polynomials p=a+b
% The polynomials are represented with the HIGHEST coefficient first
% e.g., x^2 + 3x + 4 -- [1 3 4],   x^2 + 3x --- [ 1 3 0]


k = length(b) - length(a);
p = [zeros(1,k) a] + [zeros(1,-k) b];
if(all(p==0))
  p = 0;
else
  p = p(find((p==0)==0):end);   % get rid of leading zeros
end
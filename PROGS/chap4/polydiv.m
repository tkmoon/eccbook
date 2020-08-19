function [q,r] = polydiv(a,b) 
% function [q,r] = polydiv(a,b)
% divide a(x)/b(x), and return quotient and remainder in q and r
%
% The polynomials are represented with the HIGHEST coefficient first
% e.g., x^2 + 3x + 4 -- [1 3 4],   x^2 + 3x --- [ 1 3 0]

 
m = length(a);
n = length(b);
q = 0;
if(length(a)<length(b))
  r = a;
else
  for j=1:m-n+1
    q(j) = a(j)/b(1);
    for l=2:n
      a(l+j-1) = a(l+j-1) - q(j)*b(l);
    end
  end
  r = a(m-n+2:m);
  if(all(r==0))
    r = 0;
  else
    r = r(find((r==0)==0):end);         % get rid of leading zeros
  end
end
function wtdist = reedsolwt(n,k)
% function reedsolwt(n,k)
%
% Compute the weight distribution for an (n,k) Reed-Solomon code
% Assumes, for convenience, that q=2^m for some m.

% Copyright, Todd K. Moon, November 2003

% n = 2^m-1
q = n+1;
dmin = n-k+1;
wtdist = zeros(q,1);
wtdist(1) = 1;
for j=1:n
  s = 0;
  for i=0:j-dmin
	s = s + (-1)^i*nchoosek(j-1,i)*q^(j-i-dmin);
  end
  wtdist(j+1) = nchoosek(n,j)*(q-1)*s;
end
  

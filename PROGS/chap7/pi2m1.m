function M = pi2m1(Pi,s,Lmax,k)
% function M = pi2m1(Pi,s,Lmax,k)
% The Koetter-Vardy algorithm A for mapping from a 
% reliability matrix Pi to a multiplicity matrix M
% with a total number of s interpolation points

fprintf('nargin=%d\n',nargin);
if(nargin==3)
  error('if Lmax is specified, k must be also');
end
if(nargin==4)
  lspec = 1;
else
  lspec = 0;
end
Pistar = Pi;
M = zeros(size(Pi));
Cost = 0;
for i=1:s
  [p,istar] = max(Pistar(:));
  Pistar(istar) = Pi(istar)/(M(istar)+2);
  m = M(istar)+1;
  if(lspec)								% see if we want to include this
	if( sqrt(2*(Cost+m)/(k-1)) > Lmax) 	% enough -- stop
	  break;
	end
  end
  M(istar) = m;
  Cost = Cost + M(istar);
  degbound = sqrt(2*Cost/(k-1))
end
Cost


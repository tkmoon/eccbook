function p = qf(xlist)
% 
% Compute the Q function:
%
% function p = qf(x)
%   p = 1/sqrt(2pi)int_x^infty exp(-t^2/2)dt
% Vector arguments are allowed.

% Copyright 1999 by Todd K. Moon

p = zeros(size(xlist));
i = 0;
for x = xlist
  i = i+1;
  if(x < 0)
	p(i) = 1- 0.5*erfc(-x/sqrt(2));
  else
	p(i) = 0.5*erfc(x/sqrt(2));
  end
end
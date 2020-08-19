function F = fht(f)
%  function F = fht(f)
% Compute the fast Hadamard transform of the vector f,
% where the length of f must be a power of 2
% 
% This transform actually works "in place"

% Todd K. Moon, March 24, 2004

n = length(f);
m = log2(n);  % number of stages of iteration

F = f;									% reserve space for transform 
for i=0:m-1								% count over the stages
i
  skip = 2^i; 							% amount to "skip" in indexing
skip
  j = 0;
j
  while(j < n)
	j1 = j;
j1
	for k1=0:(2^i-1) 					% number adjacent H2 steps
	  
	  tmp = F(j1+1); 					% do the H2 transform
	  F(j1+1) = F(j1+1) + F(j1+skip+1); 
	  F(j1+skip+1) = tmp - F(j1+skip+1);
	  j1 = j1+1;
	end
	j = j + 2*skip;
  end
end
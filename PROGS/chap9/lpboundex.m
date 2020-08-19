% Do the linear programming bound example

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

q = 2;									% alphabet size
n = 14;
krawtchouk;								% compute the polynomials symbolically
d = 6;
A0 = 1;
A = ones(1,n);							% distance distribution A_1..A_n
if(q==2)
  A(1,1:2:n) = 0;						% set odd weights to 0
end
A(1,1:d-1) = 0;							% set the distances to 0;
xidx = find(A ~= 0);					% A values to look for

% clear Amat, bmat
for k=1:ceil(n/2) 						% build the matrix of constraints
  i1 = 0;
  for i=xidx							% values to include in the sum
	i1 = i1+1;
	Amat(k,i1) = round(subs(K(k+1),'x',i));
	bmat(k) = round(subs(K(k+1),'x',0));
  end
end

% the simplex1 function expects _equality_ constraints.  To solve our
% problem, which is set up in terms of _inequality_ contstraints, introduce
% surplus variables.
% For example, an inequality of the form
%  4x_1 + 3x_2 >= -4,   x_1 >=0,  x_2 >= 0
% can be written as
%  4x_1 + 3x_2 - y_1 >= -4,    x_1 >= 0,  x_2 >=0,  y_1 >= 0
%
% To add these surplus variables, Amat is adjoined with a -identity,
% the c vector is adjoined to a vector of 0s.
%
% Also, since simplex1 want to _minimize_, while we want to maximize,
% we use a -c vector.
% 
% The bmat passed in is obtained, of course, by simply moving
% the column of constant terms to the other side.

Amat = [Amat -eye(size(Amat,1))];		% introduce slack variables
c = [ones(1,length(xidx)) zeros(1,size(Amat,1))]';
[x,value] = simplex1(Amat,-bmat',-c);

function [x,value,w] = simplex1(A,b,c,freevars)
% 
% Find the solution of a linear programming problem in standard form
%  minimize c'x
%  subject to Ax=b
%             x >= 0
%
% function [x,value,w] = simplex1(A,b,c,freevars)
% 
% A, b, c: system problem
% freevars = (optional) list of free variables in problem
%
% x = solution
% value = value of solution
% w = (optional) solution of the dual problem.
%    If w is used as a return value, then the dual problem is also solved.
%    (In this implementation, the dual problem cannot be solved when free
%    variables are employed.)

% Copyright 1999 by Todd K. Moon

[m,n] = size(A);
nvars = n;                              % save this in case it changes
if(m >= n)
  error('must have more variables than constraints');
end
if(rank(A) < m)
  error('degenerate matrix');
end
value = 0;                              % value of the tableau
nfree = 0;                              % number of free variables

if(nargin == 4)         % a list of free variables was passed
  [A,b,c,value,savefree,nfree] = reducefree(A,b,c,freevars);
  [m,n] = size(A);
end

% Phase I: Find a basic solution by the use of artificial variables
idx = b<0; A(idx,:) = -A(idx,:); b(idx) = -b(idx);
tableau = [A eye(m) b;  -sum(A,1) zeros(1,m) -sum(b)];
tableau
[mn,nn] = size(tableau);
basicptr = [n+1:n+m];
[tableau,basicptr] = pivottableau(tableau,basicptr);
sbasicptr = basicptr;
B1i = tableau(1:m,n+1:n+m);             % for dual
% Build the tableau for phase II
tableau = [tableau(1:m,1:n) tableau(1:m,nn); c' value];
if(any(sbasicptr)>n)
  error('Infeasible set: cannot find initial basic feasible solution');
end
ci = tableau(end,sbasicptr)';           % for dual
% transform so there are zeros in the basic columns
for i = 1:m
  tableau(mn,:) = tableau(mn,:) - c(basicptr(i))*tableau(i,:);
end

% Phase II
[tableau,basicptr] = pivottableau(tableau,basicptr);
cf = tableau(end,sbasicptr)';           % for dual
x = zeros(1,n);
x(basicptr) = tableau(1:m,end);
value = -tableau(end,end);

if(nfree)
  x = restorefree(x,savefree,freevars);
end

if(nargout==3)
  if(nargin == 4)
    error('Cannot find dual with free variables');
  end
  w = B1i'*(ci-cf);                     % fix solution
end
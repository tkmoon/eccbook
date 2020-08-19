function [A,b,c,value,savefree,nfree] = reducefree(A,b,c,freevars)
% 
% Perform elimination on the free variables in a linear programming problem
%
% function [A,b,c,value,savefree,nfree] = reducefree(A,b,c,freevars)
% 
% A,b,c = parameters from linear programming problem
% freevars = list of free variables
%
% A,b,c = new parameters for linear programming problem (with free variables
%         eliminated)
% value = value of linear program
% savefree = tableau information for restoring free variables
% nfree = number of free variables found

% Copyright 1999 by Todd K. Moon

nfree = 0;
[m,n] = size(A);
tableau = [A b; c' 0];
[mn,nn] = size(tableau);
rowout = logical(zeros(1,mn));
savefree = [];
nvars = n;
for i=1:n
  if(freevars(i))                       % pivot on this column
    nfree = nfree+1;
    [p,idx] = max(abs(tableau(1:m,i)));
    % pivot on row idx
    if(p==0)
      error('degenerate matrix')
    end
    rowout(idx) = 1;
    p = tableau(idx,i);
    tableau(idx,:) = tableau(idx,:)/p;
    savefree = [savefree; tableau(idx,:)];
    for k=1:mn
      if(~rowout(k))
        tableau(k,:) = tableau(k,:) - tableau(idx,:) .* tableau(k,i);
      end
    end
  end
end
rowout = rowout(1:end-1);
b = -tableau(~rowout,end);
A = -tableau(~rowout,~freevars);
c = tableau(mn,~freevars)';
value = tableau(mn,nn);
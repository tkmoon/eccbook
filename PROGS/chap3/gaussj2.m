function [success,Ainv,Ac,pidx,Asyst] = gaussj2(A)
% function [success,Ainv,Ac,pidx,Asyst] = gaussj2(A)
%
% Do Gaussian elimination over GF(2) on the A matrix 
%
% success = 0 if matrix is singular
% success = 1 if matrix is nonsingular
%
% Ainv is the "inverse" of A (after A has been appropriately column
%   permuted
% Ac is the parity check part of A
% Asyst is the systematic reprsentation of A (optional return argument)
%

% This is not exactly efficient, since it deals with elements as doubles.
% A more efficient way would be to pack bits into bytes, and operate
% on multiple columns simultaneously...

[m,n] = size(A);
if(m > n)
  error('Matrix should be square or wide');
end

A1 = [A eye(m)];

success = 1;
pidx = 1:n;
sing = 0;
nstep = 0;				% keep track of the number of elim steps
for i=1:m  % loop over the columns of A
  % find the first nonzero element in row i
  for k=i:n
    if(A1(i,k)) break; end;
  end
  if(A1(i,k)==0)
    success = 0;
    % error('matrix singular in gaussj2');
    fprintf(1,'matrix singular in gaussj2: i=%d\n',i);
    sing = 1;
    break;
  end
  if(sing) break; end;
  % if k ~= i, then swap columns
  if(k ~= i)
    tmp = A1(:,i);  A1(:,i) = A1(:,k);  A1(:,k) = tmp;
    t = pidx(i);  pidx(i) = pidx(k); pidx(k) = t;
  end
  % reduce down the columns
  for k=i+1:m
	if(A1(k,i)) % if nonzero in this position, do row op to cancel
	  % rowk <-- rowi + rowk
	  A1(k,i:end) = mod(A1(i,i:end) + A1(k,i:end),2);
      nstep = nstep + 1;
    end
  end
end
% now work back up
for i=m:-1:2
  for k=i-1:-1:1
    if(A1(k,i))  % if nonzero in this position, do row op to cancel
      % row k <-- rowi + rowk
      A1(k,i:end) = mod(A1(i,i:end) + A1(k,i:end),2);
      nstep = nstep + 1;
    end
  end
end
Ainv = A1(:,n+1:n+m);
Ac = A1(:,m+1:n);
if(nargout==5)
  Asyst = A1(:,1:n);
end

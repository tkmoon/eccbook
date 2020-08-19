function [c] = masseymodM(y,M)
%
% Return the shortest binary (GF(2)) LFSR consistent with the data sequence y
%
% function [c] = massey(y)
%
% y = input sequence 
%
% c = LFSR connections, c = 1 + c(2)D + c(3)D^2 + ... c(L+1)D^L
%     (Note: opposite from usual Matlab order)

% Copyright 1999 by Todd K. Moon

N = length(y);
% Initialize the variables
Ln = 0;      % current length of LFSR
Lm = 0;      % length before last change
c = 1;       % feedback connections
p = 1;       % c before last change
s = 1;       % amount of shift

for n=1:N    % N = current matching output sequence length
  d = mod(c*y(n:-1:n-Ln)',M);     % compute the discrepancy (binary arith.)
  if(d == 0)                      % no discrepancy
    s = s+1;
  else
    if(2*Ln > n-1)                % no length change in update
      c = mod(c + [zeros(1,s) p zeros(1,Ln-(Lm+s))],M);
      s = s+1;
    else                          % update with new length
      t = c;
      c = mod([c zeros(1,Lm+s-Ln)] + [zeros(1,s) p],M);
      Lm = Ln;  Ln = n - Ln;      p = t;   s = 1;
    end
  end
end
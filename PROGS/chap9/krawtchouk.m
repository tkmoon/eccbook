% Explicitly compute a set of Krawtchouk polynomials
% Todd K. Moon, Sept 2, 2004

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

n = 14;      % set these parameters when calling this function
q = 2;

Nmax = 14;   % number of polynomials to find
syms x
K = sym('K');
K(1) = 1;             % K_0
K(2) = n*(q-1)-q*x;   % K_1
for k=1:Nmax
  K(k+2)=simplify(1/(k+1)*((k + (q-1)*(n-k) - q*x)*K(k+1)-(q-1)*(n-k+1)*K(k)));
  % compute K_{k+1}
end
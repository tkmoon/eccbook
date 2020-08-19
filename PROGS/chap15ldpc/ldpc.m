% Demonstrate LDPC decoding

% A used in chapter
A = [
1 1 1 0 0 1 1 0 0 1
1 0 1 0 1 1 0 1 1 0
0 0 1 1 1 0 1 0 1 1
0 1 0 1 1 1 0 1 0 1
1 1 0 1 0 0 1 1 1 0];

% Inverse of first part, to get systematic form (not needed for decoding)
Apinv = [1 0 1 1 0;
	     0 1 1 0 1;
		 0 1 0 1 1;
		 1 1 0 1 0;
		 1 0 1 0 1];
H = mod(Apinv*A,2);   % systematic parity check matrix

[M,N] = size(A);
K = N-M;
P = H(:,N-K+1:N);
G = [P; eye(K)];  % now A*G = 0 (mod 2)

m = [1 0 1 0 1]';
c = mod(G*m,2);
t = 2*(2*c-1);

a = 2;   % signal amplitude
sigma2 = 2;  % noise variance

% First set the channel posterior probabilities
p1 = [.22  .16  .19  .48 .55  .87 .18 .79 .25 .76]; 

% then compute the received values that correspond to these
r =  log((1./p1)-1)/(-2*a)*sigma2;  % received vector

x0 = p1 > 0.5;
z0 = mod(A*x0',2);

Nloop = 50;

Lc = 2*a/sigma2;

x = probdecode(A,p1,Nloop)

fprintf('prob decoding done\n');

logrrat = -Lc*r;   % channel likelihood ratios

% xlog = ldpclogdecode(A,logrrat,Nloop,Lc)


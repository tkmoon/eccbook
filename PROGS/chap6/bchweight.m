% Get the weight distribution, starting from the weight distribution of the
% dual code

syms z A B

m = 5;
n=2^m-1;
i1 = 2^(m-1)-2^((m+1)/2-1);
n1 = (2^(m-2) + 2^((m-3)/2))*(2^m-1);
i2 = 2^(m-1);
n2 = (2^(m-1)+1)*(2^m-1);
i3 = 2^(m-1)+2^((m+1)/2-1);
n3 = (2^(m-2) - 2^((m-3)/2))*(2^m-1);
B = 1 + n1*z^i1 + n2*z^i2 + n3*z^i3;

k= n -log2(subs(B,z,1));  % weight of code described by A
% B = 1 + 310*z^12  + 527*z^16 + 186*z^20;
A = simplify(expand((1+z)^n* subs(B,z,(1-z)/(1+z))))/2^(n-k);
Avec = sym2poly(A);  % get the vector of coefficients

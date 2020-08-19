
p00 = 0.98;
p10 = 0.02;
p01 = 0.003;
p11 = 0.997;

s  = 0:0.001:0.999;
Z = (p00.^(1-s)) .* (p01.^s) + (p10.^(1-s)) .* (p11.^s);

plot(s,Z);
[Zmin,i] = min(Z);
Zmin
s(i)


syms T D N
k = 1;
bdfree =1;
dfree = 5;
T = D^5*N/(1-2*D*N);
Pej = subs(subs(T,N,1),D,Zmin);
Pb = subs(subs(diff(T,N),N,1),D,Zmin);
Pbapprox = bdfree*Zmin^dfree
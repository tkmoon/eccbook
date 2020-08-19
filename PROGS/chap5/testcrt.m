% An example of CRT calculations
m = [4 27 25];
gamma = crtgamma(m);  % precompute the gamma values

a = [0 2 3];
x = fromcrt(a,m);


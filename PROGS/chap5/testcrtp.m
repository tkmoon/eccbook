% set up data to demonstrate CRT over polynomials
m{1} = [1 -1];   m{2} = [1 -4 4];   m{3} = [1 -9 27 -27];

mp = conv(conv(m{1},m{2}),m{3});

f1 = [1 4 5 2 3 2];
a1 = tocrtpoly(f1,m);

[f1,gamma] = fromcrtpoly(a1,m);

f2 = [1 2 0 0];
a2 = tocrtpoly(f2,m);

a1a2 = tocrtpoly(polyadd(f1,f2),m);
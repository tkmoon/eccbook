
% Some generator matrices and lattice coding gains
% for the lattices A2, D4, E6, E8, Lambda16, Lambda24


A2gen = [1 0; .5 sqrt(3)/2];
d2 = 1;
v2 = sqrt(det(A2gen*A2gen'));
v2_d = v2^(1/2);
g2 = 1/v2_d^2;
tau2 = 6;   % kissing number

A3gen = [-1 -1 0; 1 -1 0; 0 1 -1];
d3 = sqrt(2);
A3gen = A3gen/d3;  % normalize so min dist is 1.
v3 = sqrt(det(A3gen*A3gen'));  % fundamental volume
v3_d = v3^(1/3);               % corresponding distance
g3 = 1/v3_d^2;                 % lattice coding gain
tau3 = 12;   % kissing number




D4gen = [1 0 0 0; 0 1 0 0; 0 0 1 0; .5 .5 .5 .5];
dmin4 = 1;
v4 = sqrt(det(D4gen*D4gen'));  % fundamental volume
v4_d = v4^(1/4);               % corresponding distance
g4 = 1/v4_d^2;				   % lattice coding gain
tau4 = 24;

E6gen = [
0 -1 1 0 0 0 0 0 
0 0 -1 1 0 0 0 0
0 0 0 -1 1 0 0 0
0 0 0 0 -1 1 0 0
0 0 0 0 0 -1 1 0
.5 .5 .5 .5 -.5 -.5 -.5 -.5];
dmin6 = sqrt(2);
E6gen = E6gen/dmin6;  % normalize so minimum distance in lattice=1
v6 = sqrt(det(E6gen*E6gen'));  % fundamental volume
v6_d = v6^(1/6);   % equivalent edge length
g6 = 1/v6_d^2;     % energy from equivalent edge length
tau6 = 72;
	
E8gen = [2 0 0 0 0 0 0 0 
	-1 1 0 0 0 0 0 0 
	0 -1 1 0 0 0 0 0  
	0 0 -1 1 0 0 0 0 
	0 0 0 -1 1 0 0 0
	0 0 0 0 -1 1 0 0
	0 0 0 0 0 -1 1 0
    .5 .5 .5 .5 .5 .5 .5 .5];
dmin8 = sqrt(2);
E8gen = E8gen/dmin8;
v8 = sqrt(det(E8gen*E8gen'));
v8_d = v8^(1/8);
g8 = 1/v8_d^2;
tau8 = 240;

L16gen = [
4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0
2 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0
2 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0
2 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0
2 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0
2 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0
2 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0
2 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0
2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0
2 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0
1 1 1 1 0 1 0 1 1 0 0 1 0 0 0 0
0 1 1 1 1 0 1 0 1 1 0 0 1 0 0 0
0 0 1 1 1 1 0 1 0 1 1 0 0 1 0 0
0 0 0 1 1 1 1 0 1 0 1 1 0 0 1 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ];  % /sqrt(2);
% The last five rows are a generator for the first-order Reed-Muller
% code of length 16.
dmin16 = sqrt(8);
L16gen = L16gen/dmin16;
v16 = sqrt(det(L16gen*L16gen'));
v16_d = v16^(1/16);
g16 = 1/v16_d^2;
tau16 = 4320;

b1 = [8 0 0 0; 4 4 0 0; 4 0 4 0; 4 0 0 4];
b2 = [4 0 0 0; 4 0 0 0; 4 0 0 0; 2 2 2 2];
b3 = [4 0 0 0; 0 4 0 0; 0 0 4 0; 2 2 2 2];
b4 = [4 0 0 0; 2 2 0 0; 2 0 2 0; 2 0 0 2];
b5 = [0 0 0 0; 2 2 0 0; 2 0 2 0; 2 0 0 2];
b6 = [4 0 0 0; 2 2 0 0; 2 0 2 0; 2 0 0 2];
b7 = [4 0 0 0; 2 0 2 0; 2 0 0 2; 2 2 0 0];
b8 = [0 0 0 0; 2 0 0 2; 2 2 0 0; 2 0 2 0];
b9 = [0 0 0 0; 2 2 0 0; 2 0 2 0; 2 0 0 2];
b10= [0 2 2 2; 0 0 0 0; 0 0 0 0; -3 1 1 1];
b11= [2 0 0 0; 0 0 0 0; 0 0 0 0; 1 1 1 1];
b12= [2 0 0 0; 2 2 0 0; 2 0 2 0; 1 1 1 1];

Z = zeros(4,4);
L24gen = [
b1 Z Z Z Z Z;
b2 b3 Z Z Z Z;
b2 Z b3 Z Z Z;
b4 b5 b5 b6 Z Z;
b7 b8 b9 Z b6 Z;
b10 b11 b12 b12 b12 b12];  % /sqrt(8);
dmin24 = sqrt(32);
L24gen = L24gen/dmin24;
v24 = sqrt(det(L24gen*L24gen'));
v24_d = v24^(1/24);
g24 = 1/v24_d^2;
tau24 = 196560;




	

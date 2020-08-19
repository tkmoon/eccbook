% Test the Feng-Tzeng algorithm

pt1 = [1,1];  pt2 = [2,3];  pt3 = [4,2];

A = [
1 pt1(1)   pt1(1)^2 pt1(1)^3   pt1(2) pt1(1)^4 pt1(1)*pt1(2) pt1(1)^5 
1 pt2(1)   pt2(1)^2 pt2(1)^3   pt2(2) pt2(1)^4 pt2(1)*pt2(2) pt2(1)^5 
0 0   0   0     1 0       pt2(1)   0
0 1   2*pt2(1) 3*pt2(1)^2 0 4*pt2(1)^3   pt2(2)   5*pt2(1)^4 
1 pt3(1)   pt3(1)^2 pt3(1)^3   pt3(2) pt3(1)^4 pt3(1)*pt3(2) pt3(1)^5 
0 0   0   0     1 0       pt3(1)   0 
0 1   2*pt3(1) 3*pt3(1)^2 0 4*pt3(1)^3   pt3(2)   5*pt3(1)^4 
];


p = 5;
A = mod(A,p);
C = fengtzeng(A,p);




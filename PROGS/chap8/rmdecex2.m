% R(2,4) Reed-Muller decoding example

% Build up the basis vectors
v1 = [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
v2 = [0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1];
v3 = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1];
v4 = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];

% Received vector
r = [0 1 1 1 0 1 0 0 0 0 0 1 0 1 0 0];  % (gives decoding failure)
r = [0 1 0 1 1 0 1 1 0 0 0 1 1 0 1 1 ]; % Wicker's example

m = [0 0 0 0 1 0 0 1 0 0 0];
c = m*[ones(1,16);v4;v3;v2;v1; v3.*v4; v2.*v4; v2.*v3; v1.*v4; v1.*v3; v1.*v2];
c = mod(c,2);

r = c;
r(9) = 1-r(9);							% flip a bit
% r = [0 1 0 1 0 1 1 0 1 1 0 1 0 1 1 0 ];


% Explicitly build up the lists of bits used in the orthogonal equations
m12list = [0 1 2 3
	       4 5 6 7
		   8 9 10 11
		   12 13 14 15];
m13list = [0 1 4 5
	       2 3 6 7
		   8 9 12 13
		   10 11 14 15];
m14list = [0 1 8 9
	       2 3 10 11
		   4 5 12 13
		   6 7 14 15];
m23list = [0 2 4 6
	       1 3 5 7
		   8 10 12 14
		   9 11 13 15];
m24list = [0 2 8 10
	       1 3 9 11
		   5 7 13 15
		   4 6 12 14];
m34list = [0 4 8 12
           1 5 9 13
	       2 6 10 14
		   3 7 11 15];

m1list = [0 1; 2 3; 4 5; 6 7; 8 9; 10 11; 12 13; 14 15];
m2list = [0 2; 1 3; 4 6; 5 7; 8 10; 9 11; 12 14; 13 15];
m3list = [0 4; 1 5; 2 6; 3 7; 8 12; 9 13; 10 14; 11 15];
m4list = [0 8; 1 9; 2 10; 3 11; 4 12; 5 13; 6 14; 7 15];

for i=1:4  % build up estimates of m12
  m12hatset(i) = mod(sum(r(m12list(i,:)+1)),2);
end
s = sum(m12hatset);
if(s > 2)
  m12hat = 1;
elseif(s < 2)
  m12hat = 0;
else
  error('Decoding failure at m12');
end

for i=1:4  % build up estimates of m13
  m13hatset(i) = mod(sum(r(m13list(i,:)+1)),2);
end
s = sum(m13hatset);
if(s > 2)
  m13hat = 1;
elseif(s < 2)
  m13hat = 0;
else
  error('Decoding failure at m13');
end

for i=1:4  % build up estimates of m14
  m14hatset(i) = mod(sum(r(m14list(i,:)+1)),2);
end
s = sum(m14hatset);
if(s > 2)
  m14hat = 1;
elseif(s < 2)
  m14hat = 0;
else
  error('Decoding failure at m14');
end

for i=1:4  % build up estimates of m23
  m23hatset(i) = mod(sum(r(m23list(i,:)+1)),2);
end
s = sum(m23hatset);
if(s > 2)
  m23hat = 1;
elseif(s < 2)
  m23hat = 0;
else
  error('Decoding failure at m23');
end

for i=1:4  % build up estimates of m24
  m24hatset(i) = mod(sum(r(m24list(i,:)+1)),2);
end
s = sum(m24hatset);
if(s > 2)
  m24hat = 1;
elseif(s < 2)
  m24hat = 0;
else
  error('Decoding failure at m24');
end

for i=1:4  % build up estimates of m34
  m34hatset(i) = mod(sum(r(m34list(i,:)+1)),2);
end
s = sum(m34hatset);
if(s > 2)
  m34hat = 1;
elseif(s < 2)
  m34hat = 0;
else
  error('Decoding failure at m34');
end

% stack up the estimates for this block
m2vechat = [m34hat m24hat m23hat m14hat m13hat m12hat];
% strip off the component decoded so far
rp = mod(r - m2vechat*[ v3.*v4; 
	             v2.*v4;
				 v2.*v3;
				 v1.*v4;
				 v1.*v3;
				 v1.*v2],2);
			 
for i=1:8  % build up estimates of m1
  m1hatset(i) = mod(sum(rp(m1list(i,:)+1)),2);
end
s = sum(m1hatset);
if(s > 4)
  m1hat = 1;
elseif(s < 4)
  m1hat = 0;
else
  error('Decoding failure at m1');
end

for i=1:8  % build up estimates of m2
  m2hatset(i) = mod(sum(rp(m2list(i,:)+1)),2);
end
s = sum(m2hatset);
if(s > 4)
  m2hat = 1;
elseif(s < 4)
  m2hat = 0;
else
  error('Decoding failure at m2');
end

for i=1:8  % build up estimates of m3
  m3hatset(i) = mod(sum(rp(m3list(i,:)+1)),2);
end
s = sum(m3hatset);
if(s > 4)
  m3hat = 1;
elseif(s < 4)
  m3hat = 0;
else
  error('Decoding failure at m3');
end

for i=1:8  % build up estimates of m4
  m4hatset(i) = mod(sum(rp(m4list(i,:)+1)),2);
end
s = sum(m4hatset);
if(s > 4)
  m4hat = 1;
elseif(s < 4)
  m4hat = 0;
else
  error('Decoding failure at m4');
end

% stack up the estimates for this block
m1vechat = [m4hat m3hat m2hat m1hat];
% strip off the component decoded so far
rpp = mod(rp - m1vechat*[v4; v3;v2;v1],2);
s = sum(rpp);
if(s > 8)
  m0hat = 1;
elseif(s < 8)
  m0hat = 0;
else
  error('Decoding failure at m0');
end

% Stack up the entire decoded vector

mhat = [m0hat m1vechat m2vechat];

  

% R(1,3) Reed-Muller decoding example

m = 3;
n = 2^m;

H = [1 1 1 1 1 1 1 1;					% Hadamard transform matrix
	 1 -1 1 -1 1 -1 1 -1;
	 1 1 -1 -1 1 1 -1 -1;
	 1 -1 -1 1 1 -1 -1 1;
	 1 1 1 1 -1 -1 -1 -1;
	 1 -1 1 -1 -1 1 -1 1;
	 1 1 -1 -1 -1 -1 1 1;
	 1 -1 -1 1 -1 1 1 -1];
% Build the basis vectors
v{1} = [0 1 0 1 0 1 0 1];
v{2} = [0 0 1 1 0 0 1 1];
v{3} = [0 0 0 0 1 1 1 1];

r = [1 0 0 1 0 0 1 0 ]; 				% received signal vector
R = 1-2*r; 								% compute (-1)^r_i
T = R*H;								% Hadamard transform
[tm,i] = max(abs(T)); 					% find maximum value
i = i-1; 								% set to 0-based indexing
ibin = dec2bin(i)-'0';					% convert to binary (in numeric form)
c = zeros(1,n);
for j=1:length(ibin)					% add up all the components
  c = c + ibin(length(ibin)-j+1)*v{j};
end
c = mod(c,2);							% binary additions
if(T(i+1) < 0)							% complement the bits
  fprintf(1,'complementing\n');
  c = 1-c;
end
c
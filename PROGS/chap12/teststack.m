% Test the stack algorithm
if 1
 r = [1 1
	1 0
    0 0 % 	1 0   % 00
	1 0
	1 1
	0 1
	0 0
	0 1]';
n = 2;
k = 1;
pc = 0.1;
nextstate = [0 1
	         2 3
			 0 1
			 2 3];
output = [0 3
	      1 2
		  3 0
		  2 1];

Delta = 10;
fanoalg(r,k,n,nextstate,output,pc,Delta);

% stackalg(r,k,n,nextstate,output,pc);

return;

end

if 0
% Example from L&C
r = [0 1 0
	 0 1 0
	 0 0 1
	 1 1 0
	 1 0 0
	 1 0 1
	 0 1 1]';

n = 3; k = 1;
pc = 0.1;

nextstate = [0 1
	         2 3
			 0 1
			 2 3];
output = [0 7
	      5 2
		  3 4
		  6 1];

% stackalg(r,k,n,nextstate,output,0.1);
Delta = 1;
fanoalg(r,k,n,nextstate,output,pc,Delta);

end

if 0
% Example from Wicker
r = [1 1
	 1 1
	 1 1
	 0 1
	 1 1
	 0 0
	 1 0]';
n = 2;
k = 1;
pc = 0.125;
nextstate = [0 1
	         2 3
			 4 5
			 6 7
			 0 1
			 2 3
			 4 5
			 6 7];
output = [0 3
	      1 2
		  2 1
		  3 0
		  3 0
		  2 1
		  1 2
		  0 3];
Delta = 10;
fanoalg(r,k,n,nextstate,output,pc,Delta);
end

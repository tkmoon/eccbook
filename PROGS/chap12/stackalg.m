function a = stackalg(r,k,n,nextstate,output,pc)
%
% r = [r(0), r(1), ... r(L-1)],  where each r(j) is a column vector of length n
% k = number of input bits for conv. code
% n = number of output bits for conv. code
% nextstate(state,input) is the next state in conv. encoder
% output(state,input) is the output sequence (binary-to-decimal converted)
% pc = crossover probability for BSC

L = size(r,2);							% number of branches

S{1}.m = 0;								% path metric
S{1}.x = [];							% input sequence
S{1}.state = 1;							% state (+1)

N = 1; 									% number of elements in stack
R = k/n;								% rate of code

stepctr = 0;
while(1)
  % Take the top path on the stack, and extend it
  state = S{1}.state;					% get the state of the best path
  M = S{1}.m;							% get the metric of best path
  j = length(S{1}.x)+1; 				% number of branch being extended to
  rj = r(:,j);
% fprintf(1,'r='); fprintf(1,'%d ',rj); fprintf(1,'\n');
  for k1=0:2^k-1 						% for each possible input
	ns = nextstate(state,k1+1);			% get the next state
	outp = output(state,k1+1);			% get the outputs
	outps = dec2bin(outp,n)-'0';		% convert to bit string
	M = S{1}.m + fanomet(rj,outps,R,n,pc); % compute fano metric
	% fprintf(1,'k1=%d  outps=',k1); fprintf(1,'%d ',outps); fprintf(1,'M=%d\n',M);
	N = N+1;							% now add this to stack
	S{N}.m = M;							% metric
	S{N}.x = [S{1}.x k1];				% input sequence
	S{N}.state = ns+1;
  end
  % Now sort by metric, largest to smallest
  mlist = zeros(1,N-1);
  for i=2:N								% pull out all the metric info.
	mlist(i-1) = S{i}.m; 				% (skip first one, to be deleted)
  end
  [msort,idx] = sort(mlist);
  
  %idx(N-1) = index of largest 
  %idx(1) = index of smallest

  % rebuild in sorted order
  for i=1:N-1
	newS{i}.m = mlist(idx(N-i));
	newS{i}.x = S{idx(N-i)+1}.x;
	newS{i}.state = S{idx(N-i)+1}.state;
  end
  S = newS;
  N = N-1; 								% to account for throwing away first
  % print it out
  stepctr = stepctr+1;
  fprintf(1,'Step %d\n',stepctr);
  for i=1:N
	fprintf(1,'%.2g & [',S{i}.m);
	fprintf(1,'%d',S{i}.x);
	fprintf(1,']\\\\\n');
  end
  if(length(S{1}.x) == L)				% if number of inputs=input length
	break;								% done!
  end
% input('press return')
end
	
  
  
  
	
	
	
	
  
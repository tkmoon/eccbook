function a = fanoalg(r,k,n,nextstate,output,pc,Delta)
%
% r = [r(0), r(1), ... r(L-1)],  where each r(j) is a column vector of length n
% k = number of input bits for conv. code
% n = number of output bits for conv. code
% nextstate(state,input) is the next state in conv. encoder
% output(state,input) is the output sequence (binary-to-decimal converted)
% pc = crossover probability for BSC

global firstvisitlist;

firstvisitlist = {};

L = size(r,2);							% number of branches

T = 0;									% threshold

inpseq = zeros(1,L); 					% input sequence to current node
stateseq = zeros(1,L+1); 				% (state sequence+1) to current node
stateseq(1) = 1;						% (initial state+1) of encoder
bestest = zeros(1,L); 					% order of try: 1=best,2=nextbest, etc.
Mlist = zeros(1,L+2); 					% list of metrics
Mlist(1) = -1e99;						% metric backward from root
Mlist(2) = 0;							% metric at root

N = 1; 									% number of elements in stack
R = k/n;								% rate of code

stepctr = 0;

metlist = zeros(1,2^k);					% place to stick metrics computed 
dobreak2 = 0;

nsteps = 0;
n1 = 0;									% index to path length

while(1)  % Loop A
%fprintf(1,'loop A\n');
  best = 1;
  while(1) 								% Loop B
%fprintf(1,'loop B\n');
  nsteps = nsteps + 1;
  fprintf(1,'\n\\vbox{\\noindent$n=%d$: ',nsteps);
fprintf(1,'$T=%.2g$ ',T);  
    rj = r(:,n1+1); 					% next input
    state = stateseq(n1+1); 			% state at current node
    % get information on all forward nodes
	M = Mlist(n1+2);
    for k1=0:2^k-1 						% for each possible input
	  ns = nextstate(state,k1+1); 		% get the next state
	  outp = output(state,k1+1); 		% get the outputs
	  outps = dec2bin(outp,n)-'0'; 		% convert to bit string
	  metlist(k1+1) = M + fanomet(rj,outps,R,n,pc); % compute fano metric
	end
	[msort,idx] = sort(metlist);		% sort in _increasing_ order
P = msort(end:-1:1);
fprintf(1,'$P=['); fprintf(1,'%.2g ',P); fprintf(1,']$ ');
fprintf(1,'$t_%d=%d$\\\\\n',n1,best-1);
	% look forward to best node, then next best, etc.
	kbest = idx(end-best+1)-1; 		% best value of k
    Mf = metlist(idx(end-best+1)); 	% largest metric (looking ahead)

fprintf(1,'Look forward: $M_F=%.2g$\\\\\n',Mf);
    if(Mf >= T) 							% Metric exceeds threshold
fprintf(1,'$M_F \\geq T$.  Move forward\\\\\n');
      n1 = n1+1;
      Mlist(n1+2) = Mf;					% save the metric
	  inpseq(n1) = kbest;				% save the input
	  stateseq(n1+1) = nextstate(state,kbest+1)+1; % save the state
	  bestest(n1) = best;				% save the branch number
	  Mb = Mlist(n1+1);					% the metric looking back (to print)
%pr(inpseq);
%pr(stateseq);
%pr(Mb);
%pr(Mlist);
%pr(bestest);

	  if(n1 == L) 			% we have reached the end of sequence
		dobreak2 = 1; 					% done!
	  end
	  if(~firstvisit(inpseq(1:n1))) 	% if not first visit, loop back
% fprintf(1,'not first visit\n');
		break; 							% loop to A
	  else 								% if first visit
fprintf(1,'First visit:  Tighten $T$\\\\\n');
		while(T+Delta <= Mf) 			% tighten the threshold
		  T = T + Delta;
		end
		break;							% loop to A
	  end
	else 								% Metric does not exceed threshold
fprintf(1,'$M_F < T$: Look back\\\\\n');	  
	  dobreak = 0;
	  while(1)							% loop C
%fprintf(1,'loop C\n');
        Mb = Mlist(n1+1);				% metric one stage back in path
fprintf(1,'$M_B=%.2g$\\\\\n',Mb);
		if(Mb <  T) 					% look back
fprintf(1,'$M_B < T$: $T = T-\\Delta$\\\\ \n');
		T = T-Delta;
		  dobreak = 1;
	      break;						% loop to A
		else 							% Mb >= T; move back
fprintf(1,'$M_B \\geq T$: Move back\\\\\n');
		  best = bestest(n1);
          n1 = n1-1;
%		  Mf = Mlist(n1+2);
%pr(Mlist);
%fprintf(1,'Mf(back)=%.2g  best(back)=%d\n',Mf,best);

		  if(best == 2^k) 				% worst node: look back
fprintf(1,'No more forward nodes\\\\\n');
			continue; 					% loop to C
		  end
		  best = best+1; 				% 1st choice, then second, etc.
fprintf(1,'All forward nodes not yet tested.  $t_%d=%d$\\\\\n',n1,best-1);
		  % else: look forward to next best node
		  break; 						% loop to B (since dobreak = 0)
		  % 
		end
	  end % while(1) loop C
	  if(dobreak) break; end;   % break out of loop B to loop A
	end
  fprintf(1,'$M_F=%.2g$ $M_B=%.2g$ Node=[',Mf,Mb);
  fprintf(1,'%d',inpseq(1:n1)); fprintf(1,']\\\\$M=%.2g$ $T=%.2g$}\\medskip\n',Mlist(n1+2),T);
	if(dobreak2) break; end;
  end % while(1) loop B
  fprintf(1,'$M_F=%.2g$ $M_B=%.2g$ Node=[',Mf,Mb);
  fprintf(1,'%d',inpseq(1:n1)); fprintf(1,']\\\\$M=%.2g$ $T=%.2g$}\\medskip\n',Mlist(n1+2),T);
  if(dobreak2) break; end;
end % while(1) loop A



function fv = firstvisit(inpseq)
% see if this is the first visit to this node.
% This is a very inefficient way to do this --- it is hard to 
% think of a way to do it worse.  But it is quick and easy 
% to code.  Perhaps a hash would be better in a real implementation
% 

global firstvisitlist;

ni = length(inpseq);

nvb = length(firstvisitlist);
for i = 1:nvb
  if(length(firstvisitlist{i}) == ni && all(firstvisitlist{i} == inpseq))
	fv = 0;
	return;
  end
end

firstvisitlist{nvb+1} = inpseq;
fv = 1;
return;

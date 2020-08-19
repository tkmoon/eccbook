function chat = ldpclogdecode(A,lambda,Nloop,Lc) 
% function chat = loggaldecode(A,r,Nloop,Lc) 
%
% Do log-likelihood decoding on a low-density parity check code
% A = parity check matrix
% lambda = channel log likelihoods = log p(c=0|y)/p(c=1|y)
% Nloop = number of iterations
% Lc = channel reliability

[M,N] = size(A);
clear Nl Ml
for m=1:M Nl{m} = []; end;
for n=1:N Ml{n} = []; end;
% Build the sparse representation of A using the M and N sets
for m=1:M
  for n=1:N
	if(A(m,n))
	  Nl{m} = [Nl{m} n];
	  Ml{n} = [Ml{n} m];
	end
  end
end

% Initialize the probabilities
Lmn = zeros(M,N);   % messages from check node m to variable node n
Lnm = zeros(N,M);   % messages from variable node n to check node m
fprintf(1,'lambda[0]:'); splatexform(1,lambda,1);

% Build the initial matrix
for n=1:N   % for each variable node
   for m = Ml{n}   % for each check node adjenct to this variable node
      Lnm(n,m) = lambda(n);
   end
end
fprintf(1,'Lnm(init)');
splatexform(1,Lnm',1);

for loop = 1:Nloop
  fprintf(1,'loop=%d\n',loop);
  for m = 1:M 			  % for each row (check)
	for n=Nl{m} % work across the columns ("horizontally")
	  pr = 1;
	  for np = Nl{m}
		if(np == n) continue; end;  % skip this for extrinsic
		pr = pr*tanh(Lnm(np,m)/2); % accumulate the product
	  end % for np
	  Lmn(m,n) = 2*atanh(pr);
	end % for n
  end % for m

  for n=1:N 			        % for each column (bit)
     Lnout(n)  = lambda(n);  % the "output" likelihood
     for m = Ml{n}	% work down the rows ("vertically")
        Lnout(n)  = Lnout(n) + Lmn(m,n);  % cumulate output likelihood
     end
     for m = Ml{n}  % remove info (extrinsic)
        Lnm(n,m) = Lnout(n) - Lmn(m,n);
     end
  end % for n

% $$$   for n=1:N 			        % for each column (bit)
% $$$      Lnout(n)  = lambda(n);  % the "output" likelihood
% $$$      for m = Ml{n}	% work down the rows ("vertically")
% $$$         Lnout(n)  = Lnout(n) + Lmn(m,n);  % cumulate output likelihood
% $$$         Lnm(n,m) = lambda(n);
% $$$         for mp = Ml{n}
% $$$            if(mp == m)  
% $$$               continue;  % otherwise, skip this for extrinsic
% $$$            end
% $$$            Lnm(n,m) = Lnm(n,m) + Lmn(mp,n);
% $$$         end % for mp
% $$$      end % for m
% $$$   end % for n
  
fprintf(1,'Lmn(%d)',loop); splatexform(1,Lmn,1);
fprintf(1,'Lnm(%d)',loop); splatexform(1,Lnm',1);
fprintf(1,'Lnout(%d)',loop); splatexform(1,Lnout,1);
%   p = 1 ./ (1+exp(lambda));  % needed only for comparison purposes!

  chat = Lnout <= 0;  % compute decoded bits for stopping criterion
  z1 = mod(A*chat',2)';

fprintf(1,'chat: ');latexform(1,chat,1);
fprintf(1,'z: ');latexform(1,z1,1);
  if(all(z1==0)) break; end
end  % end for loop

if(~all(z1==0))
  fprintf(1,'Decoding failure after %d iterations',Nloop);
end

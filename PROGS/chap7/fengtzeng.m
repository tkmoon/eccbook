function Cout = fengtzeng(A,p)
% function Cout = fengtzeng(A,p)
% Find the connection polynomial C = 1 + c1 x + c2 x^2 + ... + cl x^l
% of minimal length such that the first l+1 columns are linearly 
% dependent
% Do operations mod p (p = 2 is the default)

% Coded by Todd K. Moon, February 19, 2004

if(nargin == 1)
  p = 2;
end

demo=1;									% stuff to get info for example
% (comment out stuff after if(demo))

[M,N] = size(A);
% The McEliece representation
rho = zeros(1,N); 						% store location of nonzero discrep.
dsave = zeros(1,N);
C = 1;

s = 0;
while(1)								% loop over columns
  s = s+1;								% increment column counter
  r = 0;								% reset row counter
  columnblocked = 0;
  while(1)								% loop over rows
	r = r+1;							% next row
	drs = mod(A(r,1:s)*[zeros(s-length(C),1); C],p);  % compute discrepancy
	if(drs ~= 0) 						% nonzero discrepancy: update
	  ul = find(rho==r);			% find location with row = r
	  if(ul) 
		u = ul(1);
		l2 = length(Cs{u})+(s-u);
		l1 = length(C);
		C = [zeros(l2-l1,1); C]-(drs * invmodp(dsave(u),p))* ...
			   [zeros(l1-l2,1); Cs{u};zeros(s-u,1)];
		C = mod(C,p);
		% Now discrepancy = 0 for 1 <= i <= r
		if(demo) 
		  D(r,s) = -drs; % save, but negative to show fixed
		  fprintf(1,'Update polynomial at (%d,%d): ',r,s);
		  fprintf(1,'%d ',C); fprintf(1,'\n');
		end;
	  else  % there is nothing previous on this row
		rho(s) = r;						% save location of discrepancy
		Cs{s} = C;
		dsave(s) = drs;					% save the nonzero discrepancy
		columnblocked = 1;
		if(demo) D(r,s) = drs; end;
	  end
	end
	if(r >= M | columnblocked == 1) break; end;
  end
  if(columnblocked == 0) break; end;
end
Cout = C;
if(demo) D, rho, Cs{1:length(Cs)}, end;
	  
  

% The original Feng-Tzeng representation
% D = zeros(M,N);
% 
% s = 1;  % column counter
% r = 1;  % row counter
% C = 1;  % Stored as [c(s) c(s-1) ... c(1) c(0)=1]'
% 
% while(1)
%   drs = mod(A(r,1:s)*[zeros(s-length(C),1); C],p);  % compute discrepancy
%   if(drs==0)  % if no discrepancy
% 	if(r==M)
% 	  Cout = C;
% 	  break;
% 	else
% 	  r = r+1;
% 	end
%   else			% nonzero discrepancy
% 	% find a dru that is nonzero (on this row)
% 	ul = find(D(r,1:s-1)~=0);			% find list of locations
% 	if(ul)								% if there are any
% 	  u = ul(end);						% take the last one
% 	  dru = D(r,u);
% 	  l2 = length(Cs{u})+(s-u);
% 	  l1 = length(C);
% 	  C = [zeros(l2-l1,1); C] - (drs/dru)*[zeros(l1-l2,1); Cs{u};zeros(s-u,1)];
% 	  C = mod(C,p);
% 	  if(r==M)
% 		l = s-1;
% 		Cout = C;
% 		break;
% 	  else
% 		r = r+1;
% 	  end
% 	else
% 	  D(r,s) = drs;
% 	  Cs{s} = C;
% 	  s = s+1;
% 	  r = 1;
% 	end
%   end
% end

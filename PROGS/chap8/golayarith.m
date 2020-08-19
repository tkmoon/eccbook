% Arithmetic golay decoder

% Todd K. Moon, March 30, 2004

% test vectors
r = [0 0 0 0  0 0 1 0  0 0 0 0  1 0 0 0  0 0 0 0  0 0 0 0];
r = [0 0 1 1  1 1 0 1  0 0 0 0  0 0 1 1  0 0 0 0  1 0 1 0];
r = [1 1 1 1  0 1 1 0  1 1 1 1  0 0 1 1  0 0 0 0  1 1 1 1];


% generator matrix, systematic form.
G = [
1 0 0 0 0 0 0 0 0 0 0 0   0 1 1 1 1 1 1 1 1 1 1 1
0 1 0 0 0 0 0 0 0 0 0 0   1 1 1 0 1 1 1 0 0 0 1 0
0 0 1 0 0 0 0 0 0 0 0 0   1 1 0 1 1 1 0 0 0 1 0 1
0 0 0 1 0 0 0 0 0 0 0 0   1 0 1 1 1 0 0 0 1 0 1 1
0 0 0 0 1 0 0 0 0 0 0 0   1 1 1 1 0 0 0 1 0 1 1 0
0 0 0 0 0 1 0 0 0 0 0 0   1 1 1 0 0 0 1 0 1 1 0 1
0 0 0 0 0 0 1 0 0 0 0 0   1 1 0 0 0 1 0 1 1 0 1 1
0 0 0 0 0 0 0 1 0 0 0 0   1 0 0 0 1 0 1 1 0 1 1 1
0 0 0 0 0 0 0 0 1 0 0 0   1 0 0 1 0 1 1 0 1 1 1 0
0 0 0 0 0 0 0 0 0 1 0 0   1 0 1 0 1 1 0 1 1 1 0 0
0 0 0 0 0 0 0 0 0 0 1 0   1 1 0 1 1 0 1 1 1 0 0 0
0 0 0 0 0 0 0 0 0 0 0 1   1 0 1 1 0 1 1 1 0 0 0 1];

B = G(:,13:end);

% Now the decoder starts
s = mod(G*r',2);   % compute syndrome
if(sum(s) <= 3)							% if wt(s) <= 3
  e = [s',zeros(1,12)];
else
  found = 0;
  for i=1:12
	if(sum(mod(s + B(:,i),2)) <= 2) % if a column bi exists with wt(s+bi) <= 2
	  found = i;
	  break;
	end
  end
  if(found)
	e = [(s + B(:,i))',[zeros(1,i-1) 1 zeros(1,12-i)]];
  else
	bs = mod(B'*s,2);
	if(sum(bs) <= 3)
	  e = [zeros(12,1) bs'];
	else
	  found = 0;
	  for i=1:12
		if(sum(mod(bs + B(i,:)',2)) <= 2) %if row ri exists with wt(b's+ri)<=2
		  found = i;
		  break;
		end
	  end
	  if(found)
		e = [[zeros(1,i-1) 1 zeros(1,12-i)],bs'+B(i,:)];
	  else
		error('Too many errors');
	  end
	end
  end
end
c = mod(r + e,2);

		
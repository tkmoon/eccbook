function H = makesparseHno4(M,N,wt)
%function H = makesparseHno4(M,N,wt)
% Make a sparse parity check matrix, with column weight wt
% and row weight as uniform as possible and overlap between columns no greater
% than 1 (thus eliminating cycles of girth 4).  This is
% MacKay's construction 1A

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

H = zeros(M,N);

% Generate the matrix at random with column weight wt
for i =1:N
  colwt = 0;
  while(1)
    j = floor(M*rand)+1;
    if(H(j,i) == 0)
      H(j,i) = 1;
      colwt = colwt + 1;
    end
    if(colwt == wt) break; end;
  end
end

% Balance the row weight
rowsum = zeros(1,M);
coltoswap = 1;
while(1)
  % update row sums
  maxrs = 0; minrs = N;
  for i=1:M  % for each row
    rowsum(i) = sum(H(i,:));
    if(rowsum(i) > maxrs)
      maxrs = rowsum(i);
      maxi = i;
    end
    if(rowsum(i) < minrs)
      minrs = rowsum(i);
      mini = i;
    end
  end
  % fprintf(1,'maxrs=%d  minrs=%d\n',maxrs,minrs);
  if((maxrs - minrs) <= 1)
     break; 
  end; 					% all rows balance
  swapped = 0;
  for i=1:N   % loop to find a column to swap
    if((H(maxi,coltoswap)==1) & (H(mini,coltoswap)==0))
      H(mini,coltoswap) = 1;
      H(maxi,coltoswap) = 0;
      rowsum(maxi) = rowsum(maxi) - 1;
      rowsum(mini) = rowsum(mini) + 1;
      swapped = 1;
      break;
    end
    coltoswap = mod(coltoswap,N) + 1;
  end
  if(~swapped)
    fprintf(1,'Problem: rows not swapped:  maxrs=%d  minrs=%d\n',maxrs,minrs);
  end
end

% Now check for columns with overlap > 1
toss = [];				% list of columns to toss
for i=1:N-1
  for j=i+1:N
    overlap = sum(H(:,i) .* H(:,j));
    if(overlap > 1)  % toss either column i or j
      toss = [toss j];
    end
  end
end
toss = unique(toss);
fprintf(1,'Throwing away %d columns to satisfy no 4-loop\n',length(toss));

H(:,toss) = [];				% throw columns away
function [tableau,basicptr] = pivottableau(intableau,inbasicptr)
% 
% Perform pivoting on an augmented tableau until 
% there are no negative entries on the last row
%
% function [tableau,basicptr] = pivottableau(intableau,inbasicptr)
%
% intableau = input tableau tableau,
% inbasicptr = a list of the basic variables, such as [1 3 4]
%
% tableau = pivoted tableau 
% basicptr = new list of basic variables

% Copyright 1999 by Todd K. Moon

tableau = intableau; basicptr = inbasicptr;
[mp1,np1] = size(tableau);
n = np1-1;  m = mp1-1;

[rmin,q] = min(tableau(end,1:n));
while(rmin < 0)
  p = 0;
  minratio = realmax;
  for i=1:m
    if(tableau(i,q) > 0)
      r = tableau(i,np1)/tableau(i,q);
      if(r < minratio)
        minratio = r;
        p = i;
      end
    end
  end
  if(p == 0)
    error('unbounded solution');
  end
  % update which are the basic variables in the list
  oldb = basicptr(p); basicptr(p) = q;
  % perform the pivot
  tableau(p,:) = tableau(p,:) / tableau(p,q);
  for i=1:mp1
    if(i ~= p)
      tableau(i,:) = tableau(i,:) - tableau(p,:) .* tableau(i,q);
    end
  end
  [rmin,q] = min(tableau(end,1:n));
end
function x = restorefree(inx,savefree,freevars)
% 
% Restore the free variables by back substitution
%
% function x = restorefree(inx,savefree,freevars)
% 
% inx = linear programming solution (without free variables)
% savefree = information from tableau for substitution
% freevars = list of free variables
% 
% x = linear programming solution (including free variables)

% Copyright 1999 by Todd K. Moon

x = inx;
[nfree,nvars] = size(savefree);
nvars = nvars-1;

for i=1:nvars
  if(freevars(i))
    x = [x(1:i-1) 0 x(i:end)];
  end
end
% back substitute
j = nfree;
for i=nvars:-1:1
  if(freevars(i))
    x(i) = savefree(j,end);
    for k=1:nvars
      if(k ~= i)
        x(i) = x(i) - x(k)*savefree(j,k);
      end
    end
    j = j-1;
  end
end
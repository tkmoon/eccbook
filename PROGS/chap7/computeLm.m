function Lm = computeLm(n,k,m)
%function Lm = computeLm(n,k,m)
%
% Compute the maximum length of a list for an (n,k)
% code using GS(m) decoding

% Todd K. Moon, Feb. 12, 2004


   v = k-1;
   if(m==0) Lm = 1; return; end
   t = (v+2)/(2*v);
   Lm = floor(sqrt( n*m*(m+1)/v + t*t) - t);


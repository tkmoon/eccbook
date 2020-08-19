function x = qfinv(q)
% 
% Compute the inverse of the q function
%
% function x = qfinv(q)

% Copyright 1999 by Todd K. Moon

x = sqrt(2)*erfinv(1-2*q);


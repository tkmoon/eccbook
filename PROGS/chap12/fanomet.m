function fm = fanomet(r,oupt,R,n,pc)
% function fanomet(r,oupt,R,pc)
% 
% Compute the Fano metric for conv.-coded data through BSC
%
%r = received sequence (n symbols)
% a = coder output sequence (n symbols)
% R = rate of code
% n = number of output symbols
% pc = BSC crossover probability

r = r(:);
oupt = oupt(:); 						% stack as column vectors;
nsame = sum(r == oupt);
fm = nsame*log2(1-pc) + (n-nsame)*log2(pc) + n*(1-R);

% fm = round(fm/0.52);   % LC example
% fm = round(fm/0.31);   % Wicker example
fm = round(fm/0.348);  % my example


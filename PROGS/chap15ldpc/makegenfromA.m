function [success,G,A] = makegenfromA(A)
% function G = makegenfromA(A)
% given an MxN parity check matrix A (not necessarily systematic)
% with N>M, compute a systematic generator matrix G

% Todd K. Moon

[M,N] = size(A);
K = N-M;
[success,Ainv,Ac,pidx] = gaussj2(A);
A = A(:,pidx);
% mod(Ainv*A,2)
G = [Ac; eye(K)];

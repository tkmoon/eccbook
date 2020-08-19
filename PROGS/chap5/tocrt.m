function y = tocrt(x,m)
% function y = tocrt(x,m)
%
% Compute the representation of the scalar x using the
% using the Chinese Remainder Theorem (CRT) with
% moduli m = [m1,m2,...,mr].  It is assumed (without checking)
% that the moduli are relatively prime

y = mod(x,m);


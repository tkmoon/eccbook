
% Algebraic manipulations to derive the equations used in the
% Algebraic Golay decoder

% Uses the Matlab symbolic toolbox

syms s1 s2 s3 s4 s5 s6 s7 s8 s9 sigma1 sigma2 sigma3 
syms t1 t2 t3;
syms s32

s3 = sigma3 + sigma2*s1 + s1^3;
s32 = sigma3^2 + sigma2^2*s1^2 + s1^6;

s5 = s1^5 + sigma2*s3 + sigma3*s1^2;

s7 = s1*s32 + sigma2*s5 + sigma3*s1^4;

% Compute the numerator of D:
% use the fact that s1^3 + s3 = sigma3 + sigma2*s1
t1 =(sigma3+sigma2*s1)*(sigma3^2 + sigma2^2*s1^2) +(sigma2*s7+sigma3*s32);  
 
t1 = expand(t1);
t1 = collect(t1,sigma3);
t2 = simplify(t1/(sigma3 + sigma2*s1));
% Throw away the terms with even coefficients.
t3 = expand((sigma2 + s1^2)*(sigma2^2 + s1^4));

%Then t2 == t3.

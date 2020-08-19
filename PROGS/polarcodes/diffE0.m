

syms q0 q1 p00 p01 p10 p11  rho

%             j=0   k=0                j=0   k=1
E0 = -log(   (q0*p00^(1/(1+rho)) + q1*p01^(1/(1+rho)))^(1+rho) ...
          + (q0*p10^(1/(1+rho)) + q1*p11^(1/(1+rho)))^(1+rho) );
%             j=1  k=0          j=1   k=1
fp = diff(E0,rho);
E0diff = subs(fp,rho,0);
[n,d] = numden(E0diff);   % the denominator is equal to 1 


% k=0 j=0
I1 = q0*p00*(log(p00) - log(q0*p00 + q1*p01));

% k=0 j=1
I2 = q0*p10*(log(p10) - log(q0*p10 + q1*p11));

% k=1 j=0
I3 = q1*p01*(log(p01) - log(q0*p00 + q1*p01));

% k=1 j=1
I4 = q1*p11*(log(p11) - log(q0*p10 + q1*p11));

I = I1 + I2 + I3 + I4;

simplify(n - I)


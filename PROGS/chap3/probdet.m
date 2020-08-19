% Do an example computing probability of detection


% Information for the (31,21) BCH code
if 0
n = 31;
k = 21;
dmin = 5;
% Get the weight distribution, starting from the weight distribution of the
% dual code

syms z A B
B = 1 + 310*z^12  + 527*z^16 + 186*z^20;
A = simplify(expand((1+z)^n* subs(B,z,(1-z)/(1+z))))/2^(n-k);
Avec = sym2poly(A);  % get the vector of coefficients
end

% Information for the (15,11) Hamming code
if 1
n = 15;
k = 11;
dmin = 3;
syms z A
A = 1/(n+1)*((1+z)^n + n*(1-z)*(1-z^2)^((n-1)/2));
Avec = sym2poly(A);
end

R = k/n;
SNR = 0:10;  % SNR = 10 log10(Eb/N0)
Pu = zeros(size(SNR));
Puub = zeros(size(SNR));
Pd = zeros(size(SNR));
Pdub = zeros(size(SNR));
ps = zeros(size(SNR));
snrc = 0;
for snr = SNR
  snrc = snrc+1;
  EbN0 = 10^(snr/10);
  p = qf(sqrt(2*EbN0*R));
  ps(snrc) = 1-(1-qf(sqrt(2*EbN0)))^k;;
  for j=dmin:n
	Pu(snrc) = Pu(snrc) + Avec(j+1)*p^j*(1-p)^(n-j);
	Puub(snrc) = Puub(snrc) + nchoosek(n,j)*p^j*(1-p)^(n-j);
  end
  Pd(snrc) = 1-(1-p)^n - Pu(snrc);
  Pdub(snrc) = 1-(1-p)^n;
end
clf;
semilogy(SNR,Pu,'ro-');
hold on;
semilogy(SNR,Pd,'sg-');
semilogy(SNR,Puub,'ro--');
semilogy(SNR,Pdub,'sg--');
semilogy(SNR,ps,'hc-');
legend('P_u','P_d','P_u, upper bound','P_d, upper bound','P_u, uncoded');
xlabel('SNR(dB)');
ylabel('Probability')
input('press return to save');
print -dps probdet1.ps



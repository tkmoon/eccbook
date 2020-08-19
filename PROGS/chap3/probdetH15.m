% Do an example computing probability of detection

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
  punc = qf(sqrt(2*EbN0));
  Pdub(snrc) = 1-(1-p)^n;
end
clf;
semilogy(SNR,Pu,'ro-');
hold on;
semilogy(SNR,Pd,'sg-');
semilogy(SNR,Puub,'r+--');
semilogy(SNR,Pdub,'*g--');
semilogy(SNR,ps,'dc-');
lh = legend('P_u(E)','P_d(E)','P_u(E) upper bound','P_d(E) upper bound','P_u(E) uncoded');
xh = xlabel('E_b/N_0(dB)');
yh = ylabel('Probability');
set(gca,'fontsize',15);
set(lh,'fontsize',15);
set(xh,'fontsize',15);
set(yh,'fontsize',15);
input('press return to save');
print -dps probdet1.ps


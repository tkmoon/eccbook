% Plot probability of error for BPSK

% P0 = 0.1;  % =P(+Eb)
% P1 = 1-P0;  % P(-Eb)
startsnrdb = 0;
snrstepdb = 0.1;
endsnrdb = 10;
N0 = 2;  % sigma^2 = N0/2 = 1
sigma2 = N0/2;
sigma = sqrt(sigma2);
plist = [];
snrlist = [];
for snrdb=startsnrdb:snrstepdb:endsnrdb
  snrlist = [snrlist snrdb];
  snr = 10^(snrdb/10);   % snr not in dB
  tau = sigma2/2*log(P1/P0);
  
  Eb = N0*snr*R;
  Ebsqr = sqrt(Eb);
  
  p = qf((tau + Ebsqr)/sigma)*P1 + qf((Ebsqr - tau)/sigma)*P0;
  plist = [plist p];
end
% semilogy(snrlist,plist,'r');
  
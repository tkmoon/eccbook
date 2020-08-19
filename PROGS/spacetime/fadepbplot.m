% Plot the performance of BPSK in a fading channel

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

SNRdB = 0:30;
EbN0 = 10.^(SNRdB/10);
PBPSK = qf(sqrt(2*EbN0));
clf;
markers = 'ox+*sdv^ph';
% semilogy(SNRdB,PBPSK);
Pm1 = 0.5*(1-sqrt(EbN0 ./ (1 + EbN0)));
% subplot(2,2,1);
h = semilogy(SNRdB,Pm1,[markers(1),'-']);
set(gca,'fontsize',15);
hold on; 
mlist = [2 3 4 6 10];
P{1} = Pm1;
legstr{1} = 'm=1';
i = 2;
for m=mlist
  P{m} = 0;
  for k=0:m-1
	P{m} = P{m} + nchoosek(m-1+k,k)*(1-Pm1).^k;
  end
  P{m} = P{m} .* (Pm1.^m);
  h(i) = semilogy(SNRdB,P{m},[markers(i),'-']);
  legstr{i} = sprintf('m=%d',m);
  i = i+1;
end
axis([0 30 10^(-8) 10^(-1)]);
xh = xlabel('average \gamma_b (dB)');
ylabel('average probability of bit error');
legend(h,legstr);

% now plot the approximation
Pm1 = 1./ (4*EbN0);
semilogy(SNRdB,Pm1,[markers(1), '--']);
i = 2;
for m=mlist
  Pa = nchoosek(2*m-1,m) * Pm1.^ m;
  semilogy(SNRdB,Pa,[markers(i),'--']);
  i = i+1;
end

% print -dps fadeplot2.ps
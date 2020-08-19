
% Compute probability of error for (n,1) repetition codes, n odd

startsnrdb = 0;
snrstepdb = 0.5;
endsnrdb = 10;

markerlist = 'ox+*sd';

% plot uncoded
snrlistdb = startsnrdb:snrstepdb:endsnrdb;  
snrlist = 10 .^(snrlistdb./10);  % Eb/N0
Pbunc = qf(sqrt(2*snrlist));
figure(1);
clf
llist = semilogy(snrlistdb,Pbunc);
grid on ;
hold on;

nlist = [3, 11];
ni = 0;
for n = nlist
  ni = ni+1;
  R = 1/n;
  dmin = n;
  t = (n-1)/2;

  N0 = 2; 								% sigma^2 = N0/2 = 1
  sigma2 = N0/2;
  sigma = sqrt(sigma2);
  plist = [];
  for snrdb=startsnrdb:snrstepdb:endsnrdb
	snr = 10^(snrdb/10); 				% snr not in dB
	
	Ec = 1;
	d = 2*Ec;
	sigma2 = 1/(2*R*snr);
	sigma = sqrt(sigma2);
	p = qf(d/(2*sigma));
	
	Pe = 0;
	for i=t+1:n
	  Pe = Pe + nchoosek(n,i)*p^i*(1-p)^(n-i);
	end
	plist = [plist Pe];
  end
  
  
  semilogy(snrlistdb,plist);
  llist = [llist semilogy(snrlistdb,plist,markerlist(ni))];
  drawnow
end
xh=xlabel('E_b/N_0 (dB)');
yh = ylabel('Probability of bit error');
lh = legend(llist, 'Uncoded','Repetition n=3','Repetition n=11');

set(gca,'fontsize',15)
set(xh,'fontsize',15)
set(yh,'fontsize',15)
set(lh,'fontsize',15)

print -dps repcodeps.ps
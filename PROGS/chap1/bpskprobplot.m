% set up stuff for bpskprob

P0 = 0.5;  % =P(+Eb)
P1 = 1-P0;  % P(-Eb)
R = 1;
bpskprob;

figure(1);
clf
semilogy(snrlist,plist);
xh=xlabel('E_b/N_0 (dB)');
yh = ylabel('Probability of bit error');
grid on
set(gca,'fontsize',15)
set(xh,'fontsize',15)
set(yh,'fontsize',15)
print -dps bpskprob.ps

figure(2);
clf

semilogy(snrlist,plist);  % plot the old rate 1 stuff
grid on
drawnow 

hold on;

P0 = 0.5;  % =P(+Eb)
P1 = 1-P0;  % P(-Eb)
R = 1/2;
bpskprob;
semilogy(snrlist,plist,'--');

P0 = 0.5;  % =P(+Eb)
P1 = 1-P0;  % P(-Eb)
R = 1/3;
bpskprob;
semilogy(snrlist,plist,'-.');


xh=xlabel('E_b/N_0 (dB)');
yh = ylabel('Probability of bit error');
set(gca,'fontsize',15)
set(xh,'fontsize',15)
set(yh,'fontsize',15)
legend('Rate R=1','Rate R=1/2', 'Rate R=1/3')
print -dps bpskprobR.ps
return;

 


hold on;
P0 = 0.25;  % =P(+Eb)
P1 = 1-P0;  % P(-Eb)
bpskprob;
semilogy(snrlist,plist,'--');

hold on;
P0 = 0.1;  % =P(+Eb)
P1 = 1-P0;  % P(-Eb)
bpskprob;
semilogy(snrlist,plist,':');



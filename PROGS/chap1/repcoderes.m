
% plot results from repcode simulation program testrepcode.cc

% soft-decision decoding results
repcoderes3 = [
 0 0.0821018
0.5 0.0573394
1 0.0541712
1.5 0.0455996
2 0.0346021
2.5 0.0331455
3 0.0207082
3.5 0.0158529
4 0.0139179
4.5 0.00802439
5 0.00638081
5.5 0.00339006
6 0.00192105
6.5 0.00136435
7 0.000829469
7.5 0.000424356
8 0.000183465
8.5 7.878e-05
9 3.37826e-05
9.5 1.24943e-05
10 4.33268e-06
];


repcoderes11 = [
0 0.0653595
0.5 0.0628931
1 0.0584795
1.5 0.0505306
2 0.0391083
2.5 0.0259605
3 0.0215936
3.5 0.0175994
4 0.0130378
4.5 0.00910581
5 0.00593754
5.5 0.00369549
6 0.00265308
6.5 0.00146933
7 0.000774251
7.5 0.000478133
8 0.000201855
8.5 8.07666e-05
9 2.95929e-05
9.5 1.23578e-05
10 3.46054e-06
];


% hard-decision decoding results
repcoderesH3 = [
0 0.118483
0.5 0.0958773
1 0.0789889
1.5 0.0753012
2 0.0568828
2.5 0.0507357
3 0.0389864
3.5 0.0324675
4 0.0275558
4.5 0.0205423
5 0.0192086
5.5 0.012663
6 0.00827472
6.5 0.00579475
7 0.0031409
7.5 0.00203558
8 0.000973255
8.5 0.000644268
9 0.000335762
9.5 0.00014021
10 7.14826e-05
];


repcoderesH11 = [
0 0.121212
0.5 0.108108
1 0.0831947
1.5 0.0917431
2 0.0775194
2.5 0.0610501
3 0.0574053
3.5 0.044603
4 0.0329815
4.5 0.02886
5 0.0221631
5.5 0.0172236
6 0.0108319
6.5 0.0074096
7 0.00468867
7.5 0.00333411
8 0.00196252
8.5 0.00123267
9 0.000586961
9.5 0.000331532
10 0.000170123
];


snrlist = repcoderes3(:,1);  
snr = 10.^(snrlist/10);
plist = qf(sqrt(2*snr(:)'));

figure(1);
clf;
hlist = semilogy(snrlist,plist);  % uncoded
grid on
hold on;
% soft decoding
hlist = [hlist;semilogy(snrlist,repcoderes3(:,2),'o-')];  % rep n=3
hlist = [hlist;semilogy(snrlist,repcoderes11(:,2),'s-')];  % rep n=11
% hard decoding
hlist = [hlist;semilogy(snrlist,repcoderesH3(:,2),'o--')];  % rep n=3
hlist = [hlist;semilogy(snrlist,repcoderesH11(:,2),'s--')];  % rep n=11
hl = legend(hlist,'Uncoded','Rep n=3, soft', 'Rep n=11, soft',...
	'Rep n=3, hard','Rep n=11, hard');
xh = xlabel('SNR, dB')
yh = ylabel('Probability of bit error');
set(gca,'fontsize',15);
set(hl,'fontsize',15);
set(xh,'fontsize',15);
set(yh,'fontsize',15);
print -dps repcodeprob.ps

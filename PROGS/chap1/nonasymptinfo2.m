% Make plots for nonasymptotic information theory
%  TKM 1/3/19
% For Gaussian

nlist = [10:10:90 100:50:2000];
% nlist = [200:50:2000];
epsilon = 1e-3;
delta = 0.11;

P = 1;

C = 0.5*log2(1 + P);

V = P*(P+2)*(log2(exp(1)))^2/(2*(P+1)^2);

Qinveps = qfinv(epsilon);

MlistGauss = [];

for n = nlist
   M = C - sqrt(V/n)*Qinveps + 0.5*log2(n)/n;
   if(M < 0), M = 0; end;
   MlistGauss = [MlistGauss M];
end


figure(1);
clf
plot([nlist(1), nlist(end)],[C C],'k','linewidth',2);
hold on;
plot(nlist,MlistGauss,'g','linewidth',2);
xlabel('Block length n');
ylabel('Rate, bits/channel use');
axis([0 nlist(end) 0 0.55]);
set(gca,'fontsize',15);
h = legend('Capacity','Gaussian approximation','location','southeast');
set(h,'fontsize',15);
grid on
printcrop('nonasympgaussP1.pdf');
print -depsc nonasympgaussP1.eps



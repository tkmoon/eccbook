function z = loghist(data,itno,plotno);

bins = -50:0.25:10;

subplot(2,2,plotno);
%N = hist(log((data+realmin)./(1-data+realmin)),bins);
N = hist(log((data)./(1-data)),bins);
s = sum(N);
N = N/s;
plot(bins,N);
xl=    xlabel('Log Likelihood Ratio L');
yl=    ylabel('Probability p_L');
tl=  title(sprintf('Iteration %d',itno-1));
set(xl,'FontSize',15');
set(yl,'FontSize',15');
set(tl,'FontSize',15');


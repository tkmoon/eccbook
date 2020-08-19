% Make a plot of fading signal amplitudes

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

clf;
subplot(2,2,1);
N = 34;
g21 = jakes(N);
g22 = jakes(N);
g23 = jakes(N);
plot(10*log10(g21));
hold on;
plot(10*log10(g22),'r-.');
% plot(10*log10(g22),'-.');
xlabel('time t/T');
ylabel('|g(t)|, dB');
axis([0 200 -30 10]);
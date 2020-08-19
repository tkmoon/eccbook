function g2 = jakes(N)
% function g2 = jakes(N)
% use Jakes method to simulate the amplitude of a fading channel
% returns |g|^2, which is the SQUARE magnitude of the channel.
% 
% e.g., to plot in dB, use plot(10*log10(g2));

% Jakes method is described in _Microwave_Mobile_Communication,
% IEEE Press, 1993.

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

M = floor(0.5*(0.5*N-1));
nlist = 1:N;
mlist = 1:M;
%betan = pi*mlist/M;
betan = 2*pi*rand(1,M);
alpha = 0;
alpha = 2*pi*rand(1);
fmT = 0.1;
fnT = fmT*cos(2*pi*mlist/N);
tT = 0:200;
j = sqrt(-1);
g = sqrt(2)*(cos(alpha)*cos(2*pi*fmT*tT) + j*sin(alpha)*cos(2*pi*fmT*tT));
gi2 = M + cos(alpha)^2;
gq2 = M + sin(alpha)^2;
giq2 = sin(alpha)*cos(alpha);
for n=1:M
 g=g+2*(cos(betan(n))*cos(2*pi*fnT(n)*tT)+j*sin(betan(n))*cos(2*pi*fnT(n)*tT));
 gi2 = gi2 + cos(2*betan(n));
 gq2 = gq2 - cos(2*betan(n));
 giq2 = giq2 + 2*sin(betan(n))*cos(betan(n));
end
g2 = g.*conj(g);
g2 = g2/(gi2 + gq2);					% normalize for unit MSE

plot(10*log10(g2));

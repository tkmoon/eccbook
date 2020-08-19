% Plot the Hamming and Gilbert bounds for binary codes

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only

delta = 0:0.01:1;
delta2 = 0:0.01:0.5;
singleton = 1-delta;
hamming = 1 - h2(delta/2);
elias = 1-h2(0.5 - sqrt(0.5*(0.5-delta2)));
plotkin = 1-2*delta2;
mceliece1 = h2(0.5-sqrt(delta2 .* (1-delta2)));
% compute the McEliece2 bound
mceliece2 = zeros(1,length(delta2));
i = 0;
for d = delta2
  i = i+1;
  u = 0:.0001:1-2*d;
  al = min(1 + mcelieceg(u.^2) - mcelieceg(u.^2 + 2*d*u + 2*d));
  mceliece2(i) = al;
end

% find the smallest of the upper bounds
lowest = zeros(1,length(delta2));
i = 0;
for d = delta2
  i = i+1;
  al = min([elias(i) mceliece1(i) mceliece2(i)]);
  lowest(i) = al;
end

clf;
subplot(2,2,1);
C = [0.9 0.9 0.9];
gilbert = 1-h2(delta2);
patch([delta2, fliplr(delta2(1:end-1))],[gilbert fliplr(lowest(1:end-1))],C);
hold on
plot(delta,singleton);
hold on;
plot(delta,hamming);
plot(delta2,gilbert);
plot(delta2,plotkin)
plot(delta2,elias)
plot(delta2,mceliece1)
plot(delta2,mceliece2)

xlabel('\delta');
ylabel('\alpha(\delta)');


printcrop('bounds.pdf');
print -deps bounds.eps

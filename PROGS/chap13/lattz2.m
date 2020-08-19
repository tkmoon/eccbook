
% make a plot of the z2 lattice to discuss shaping and lattice gain

M = eye(2);

i1range = 10;
i2range = 10;
xrange = [-10,10];
yrange = [-10,10];
lattpts = [];
enlist = [];
ijlist = [];
for i1=-i1range:i1range					% walk over a range of values
  for i2 = -i2range:i2range
	lattpt = [i1 i2]*M; 				% get the lattice point
	if(lattpt(1) < xrange(1) | lattpt(1) > xrange(2) | ...
	  lattpt(2) < yrange(1) | lattpt(2) > yrange(2))
      continue;  % throw away pts outside of plot range
	end;
	lattpts = [lattpts; lattpt];
	enlist = [enlist; lattpt(1)^2 + lattpt(2)^2];
	ijlist = [ijlist; [i1 i2]];
  end
end
figure(1); clf
plot(lattpts(:,1),lattpts(:,2),'k.');
hold on;
[enlist,idx] = sort(enlist);
lattpts = lattpts(idx,:);
ijlist = ijlist(idx,:);

r16z2 = sqrt(enlist(16));
En16z2 = sum(enlist(1:16))/16;

r32z2 = sqrt(enlist(32));
En32z2 = sum(enlist(1:32))/32;

r64z2 = sqrt(enlist(64));
En64z2 = sum(enlist(1:64))/64;

r128z2 = sqrt(enlist(128));
En128z2 = sum(enlist(1:128))/128;

r256z2 = sqrt(enlist(256));
En256z2 = sum(enlist(1:256))/256;


% x = plotcirc(0,0,r16);
% plot(x(1,:),x(2,:),'--');

x = plotcirc(0,0,r32);
plot(x(1,:),x(2,:),'k--');

x = plotcirc(0,0,r64);
plot(x(1,:),x(2,:),'k--');

x = plotcirc(0,0,r128);
plot(x(1,:),x(2,:),'k--');

x = plotcirc(0,0,r256);
plot(x(1,:),x(2,:),'k--');

axis 'equal';
axis([xrange yrange]);

set(gca,'fontsize',18);

printcrop('z2constcirc.pdf');
print -deps z2constcirc.ps

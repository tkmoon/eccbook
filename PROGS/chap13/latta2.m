
% make a plot of the A2 lattice to discuss shaping and lattice gain

M = [1 0; .5 sqrt(3)/2];


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
figure(1);  clf;
plot(lattpts(:,1),lattpts(:,2),'k.');
[enlist,idx] = sort(enlist);
lattpts = lattpts(idx,:);
ijlist = ijlist(idx,:);

r16 = sqrt(enlist(16));
En16 = sum(enlist(1:16))/16;

r32 = sqrt(enlist(32));
En32 = sum(enlist(1:32))/32;

r64 = sqrt(enlist(64));
En64 = sum(enlist(1:64))/64;

r128 = sqrt(enlist(128));
En128 = sum(enlist(1:128))/128;

r256 = sqrt(enlist(256));
En256 = sum(enlist(1:256))/256;


% Now do the coset stuff
M2 = 2*M;  % generator for sublattice
lattpts1 = [];
lattpts2 = [];
lattpts3 = [];
lattpts4 = [];
for i1=-i1range:i1range					% walk over a range of values
  for i2 = -i2range:i2range
	lattpt1 = [i1 i2]*M2; 				% get the sublattice point
	if(~(lattpt1(1) < xrange(1) | lattpt1(1) > xrange(2) | ...
	  lattpt1(2) < yrange(1) | lattpt1(2) > yrange(2)))
      lattpts1 = [lattpts1; lattpt1];
	end
	lattpt2 = lattpt1 + [1 0];			% first coset point
	if(~(lattpt2(1) < xrange(1) | lattpt2(1) > xrange(2) | ...
	  lattpt2(2) < yrange(1) | lattpt2(2) > yrange(2)))
      lattpts2 = [lattpts2; lattpt2];
	end
	lattpt3 = lattpt1 + M(2,:);			% second coset point
	if(~(lattpt3(1) < xrange(1) | lattpt3(1) > xrange(2) | ...
	  lattpt3(2) < yrange(1) | lattpt3(2) > yrange(2)))
      lattpts3 = [lattpts3; lattpt3];
	end
	lattpt4 = lattpt1 + sum(M);			% third coset point
	if(~(lattpt4(1) < xrange(1) | lattpt4(1) > xrange(2) | ...
	  lattpt4(2) < yrange(1) | lattpt4(2) > yrange(2)))
      lattpts4 = [lattpts4; lattpt4];
	end
  end
end
figure(2); clf;
plot(lattpts1(:,1),lattpts1(:,2),'k.');
hold on
plot(lattpts2(:,1),lattpts2(:,2),'ko');
plot(lattpts3(:,1),lattpts3(:,2),'kx');
plot(lattpts4(:,1),lattpts4(:,2),'k+');

x = plotcirc(0,0,r16);
plot(x(1,:),x(2,:),'k--');

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

h = legend('$2\Lambda$','$2\Lambda+(1,0))$','$2\Lambda+(1/2,\sqrt{3}/2)$',...
	'$2\Lambda+(3/2,\sqrt{3}/2)$')
set(h,'interpreter','latex')
set(gca,'fontsize',18);

printcrop('a2constcirc.pdf');
print -deps a2constcirc.ps

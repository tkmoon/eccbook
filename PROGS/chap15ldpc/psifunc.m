% Make a plot of the psi function used in density evolution

% Todd K. Moon, May 6, 2004

if(~exist('psi2'))
   mlist = 1:.25:10;   %only need positive values: psi is odd
   psi = [];

   for m=mlist;
	 s2 = 2*m;
   str=sprintf('tanh(x/2).*exp(-(x-(%g)).^2 ./(4*(%g)))/sqrt(4*pi*(%g))',...
	m,m,m);
     F = inline(str);
	 Q = quad(F,m-10*s2,m+10*s2)
	 psi = [psi Q];
   end
   mlist2 = [-fliplr(mlist) 0 mlist];
   psi2 = [-fliplr(psi) 0 psi];
end

clf;
plot(mlist2,psi2);
hold on;
plot(mlist2,tanh(mlist2/2),'r--');
h = legend('\Psi(x)','tanh(x)');
set(gca,'fontsize',15)
set(h,'fontsize',15);
axis([-10 10 -1.5 1.5])
xlabel('x');
% print -dps psifunc.ps
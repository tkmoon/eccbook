
% test the algorithm for pi2m

q = 5; n = 4;
k = 2;

Lmax = 400;
Pi = rand(q,n);
for i=1:n
  Pi(:,i) = Pi(:,i)/sum(Pi(:,i));
end
Pi
s = 2000;
% M = pi2m1(Pi,s,k)
M = pi2m1(Pi,s,Lmax,k)

Mnorm = M;
for i=1:n
  Mnorm(:,i) = Mnorm(:,i)/sum(Mnorm(:,i));
end

mdat = [];
scorelist = [];
costlist = [];
slist = 10:5:400;
for s=slist
   %  M = pi2m1(Pi,s,k);
  M = pi2m1(Pi,s,Lmax,k);

  Mnorm = M;
  for i=1:n
	Mnorm(:,i) = Mnorm(:,i)/sum(Mnorm(:,i));
  end
  mdat = [mdat norm(Mnorm-Pi)];
  scorelist = [scorelist norm(M.*Pi,1)];
  costlist = [costlist sum(sum( M.*(M+1)/2,1))];
end
subplot(3,1,1);
plot(slist,mdat,'r');
ylabel('Norm Difference');
subplot(3,1,2);
plot(slist,scorelist,'g');
ylabel('Score');
subplot(3,1,3);
plot(slist,costlist,'b');
ylabel('Cost');
xlabel('Number of steps');

printcrop('mpiplot.pdf');
print -deps 'mpiplot.ps'

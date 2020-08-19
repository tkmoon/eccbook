

% Test a bunch of formulas for GS stuff

if 0
% test the formula for A(v,K) to see if it is a continuous increasing
% function of K
n = 32;
k = 8;

v = 3;
K = 0:.1:7;
r = mod(K,v);
A = K.^2/(2*v) + K/2 + r .* (v-r)/(2*v);
X = r .* (v-r)/(2*v);
plot(K,A)
return
end




if 0
% 
% Compute Lm and its bound
n = 32;
k=8;
v = k-1;
m = 120;
Lm = floor(sqrt( n*m*(m+1)/v + ((v+2)/(2*v))^2) - (v+2)/(2*v))
Lmbound = (m+0.5)*sqrt(n/v)
return
end  



if 0
% compute B
B = v*l^2/2 + (v+2)*l/2
return

end

if 1 
% Make some plots of t_m function
n=32;
k=8;
mlist = 0:125;
Kmlist = [];
for m=mlist
  Kmlist = [Kmlist computekm(n,k,m)];
end
subplot(2,1,1);
plot(mlist,Kmlist);
axis([0,max(mlist),min(Kmlist-1),max(Kmlist)+1]);
xlabel('m');
ylabel('K_m');
title(sprintf('Values of K_m for n=%d, k=%d',n,k));
printcrop('kmplot1.pdf');
print -deps kmplot1.eps
return
end


if 0

% do some computations to compute k_m and t_m

n=32;
k = 8;
fprintf(1,'m=0 t=%d\n',floor((n-k)/2));
for m=1:120
  v = k-1;
  C = n*nchoosek(m+1,2);
  k1 = 0;
  while(1)
	k1 = k1+1;
	l = m*k1;
	
	r = mod(l,v);
	n2 = l^2/(2*v) + l/2 + r*(v-r)/(2*v);

	l = l-1;
	lk = floor(l/v);
	n1 = (l+1)*(lk+1) - v/2*lk*(lk+1);
	
%	n1 = n2;
	
	if(n1 > C) break; end
  end
  K = k1;
  t = n-K;
  fprintf(1,'m=%d  K=%d  t=%d\n',m,K,t);
end
	  
end	
  
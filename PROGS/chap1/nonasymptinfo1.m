% Make plots for nonasymptotic information theory
%  TKM 1/3/19

nlist = [10:10:90 100:50:2000];
% nlist = [200:50:2000];
epsilon = 1e-3;
delta = 0.11;

P = 1;

C = 1 - H(delta);

V = delta*(1-delta)*(log2((1-delta)/delta))^2;
Qinveps = qfinv(epsilon);

MlistBSC = [];
Mlist320 = [];
Mlist310 = [];

for n = nlist
   M = C - sqrt(V/n)*Qinveps + 0.5*log2(n)/n;
   if(M < 0), M = 0; end;
   MlistBSC = [MlistBSC M];



   % Theorem 40 bound (upper bound)
   warning('OFF','MATLAB:nchoosek:LargeCoefficient');

   alpha = 1-epsilon;
   found = 0;
   for l=1:n
      alphal = 0;
      for k=0:l-1
         b = lognchoosek(n,k) + (n-k)*log(1-delta) + k*log(delta);
         alphal = alphal + exp(b);
      end
      alphas(l) = alphal;
      if(l > 1)
         if(alphas(l-1) <= alpha && alpha < alphas(l));
            L = l-1;
            alphaL = alphas(l-1);
            alphaLp1 = alphas(l);
            found = 1;
            break;
         end
      end
   end % for l
   if(found)
      lambda = (alpha - alphaL)/(alphaLp1 - alphaL);
      %       fprintf('n=%d  L=%d  lambda=%g   alpha=%g  alphaL=%g  alphaL+1=%g\n',n,L,lambda,alpha,alphaL,alphaLp1);
      betaL = 0;
      for k=0:L   % comput betaL
         b = lognchoosek(n,k) - n*log(2);
         betaL = betaL + exp(b);
      end
      b = lognchoosek(n,L+1) - n*log(2);
      betaLp1 = betaL + exp(b);
      beta = (1-lambda)*betaL + lambda*betaLp1;
      % fprintf('betaL=%g  betaLp1=%g  beta=%g\n',betaL,betaLp1,beta);
      Mthm40 = 1/beta;
      Mlist320 = [Mlist320 log2(Mthm40)/n];
   else
      error('Working lambda and L not found');
   end

   % corollary 39 bound (lower bound)

   % start the search with the upper bound found using thm. 40
   Mhi = floor(Mthm40 + .5);
   Mlo = 2;
   Mhirate = log2(Mhi)/n;
   Mlorate = log2(Mlo)/n;

   % compute all the binomial factors
   lognchooset = zeros(1,n+1);
   nchooset = zeros(1,n+1);
   for t=0:n
      lognchooset(t +1) = lognchoosek(n,t);
      nchooset(t +1) = exp(lognchooset(t +1));
   end

   sthi = computeMbound2(n,Mhi,delta,lognchooset,nchooset);   % should be < eps
   stlo = computeMbound2(n,Mlo,delta,lognchooset,nchooset);   % should be > eps

   fprintf('n=%d  Mlo=%g (%g)  stlo=%g   Mhi=%g (%g) sthi=%g\n',n,Mlo,...
                 Mlorate, stlo,Mhi,Mhirate,sthi);
   while(sthi < epsilon)
      fprintf('fixing Mhi\n');
      Mhi = 2*Mhi;
      Mhirate = log2(Mhi)/n;
      sthi = computeMbound2(n,Mhi,delta,lognchooset,nchooset);
      fprintf('n=%d  Mlo=%g (%g)  stlo=%g   Mhi=%g (%g) sthi=%g\n',n,Mlo,...
                 Mlorate, stlo,Mhi,Mhirate,sthi);
   end
   

   fprintf('---\n');

      

   %   if( (Mhi > Mlo) && sthi < epsilon && stlo > epsilon)
   if( (Mhi > Mlo) && sthi > epsilon)
      while(1)
         % search using a binary search
         Mmidrate = (Mhirate + Mlorate)/2;
         Mmid = 2^(n*Mmidrate);
         fprintf('n=%d  Mlo=%g (%g)  stlo=%g   Mhi=%g (%g) sthi=%g\n',n,Mlo,...
                 Mlorate, stlo,Mhi,Mhirate,sthi);

         stmid = computeMbound2(n,Mmid,delta,lognchooset,nchooset);
         fprintf('    Midrate=%g  stmid=%g\n',Mmidrate,stmid);
         lastMhirate = Mhirate;
         if(stmid < epsilon)
            Mlo = Mmid;
            Mlorate = Mmidrate;
            stlo = stmid;
         else
            Mhi = Mmid;
            Mhirate = Mmidrate;
            sthi = stmid;
         end
         if(stlo >= epsilon)
            Mrate = Mlorate;
            break;
         end
         if(abs(Mhirate - Mlorate) < 0.01)
            Mrate = Mlorate;
            break;
         end;
         


% $$$          if( Mhirate - Mlorate < 0.01)
% $$$             %         if(abs(Mhi - Mlo) < .1)
% $$$             M = Mhi;
% $$$             Mrate = Mhirate;
% $$$             fprintf('*** n=%d  M=%g  Mrate=%g\n',n,m,Mrate);
% $$$             break;
% $$$          end;
      end % while
   else
      Mrate = 0;
   end
   if(Mrate < 0), Mrate = 0; end;
   Mlist310 = [Mlist310 Mrate];
   
end


figure(1);
clf
plot([nlist(1), nlist(end)],[C C],'k','linewidth',2);
hold on;
plot(nlist,Mlist320,'g','linewidth',2);
plot(nlist,Mlist310,'b','linewidth',2);
plot(nlist,MlistBSC,'r','linewidth',2);
xlabel('Block length n');
ylabel('Rate, bits/channel use');
set(gca,'fontsize',15);
h = legend('Capacity','Upper bound','Lower bound',...
           'Gaussian approximation','location','southeast');
set(h,'fontsize',15);
grid on
printcrop('nonasympbsceps001.pdf');
print -depsc nonasympbsceps001.eps

function st = computeMbound(n,M,delta)
   st = 0;
   for t=0:n
      st1 = 0;
      for s=0:t
         b = lognchoosek(n,s) - n*log(2);
         st1 = st1 + exp(b);
      end
      m = min([1, (M-1)*st1]);
      b = lognchoosek(n,t) + t*log(delta) + (n-t)*log(1-delta);
      st = st + exp(b)*m;
   end % for t
end

function st = computeMbound2(n,M,delta,lognchooset, nchooset)
   st = 0;
   for t=0:n
      st1 = 0;
      for s=0:t
         b = lognchooset(s + 1) - n*log(2);
         st1 = st1 + exp(b);
      end
      m = min([1, (M-1)*st1]);
      b = lognchooset(t +1) + t*log(delta) + (n-t)*log(1-delta);
      st = st + exp(b)*m;
   end % for t
end

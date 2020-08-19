
% simple explorations of importance sampling (IS)
%  tkm 6/3/19
% d-dimensional signal constellation, interior point
%              |
%              X   
%              |
%     ----X----X----X----
%              |
%              X
%              |
% 

clear gamma
sigma = sqrt(.05);
sigma = sqrt(.02);
sigma = sqrt(.02);
sigma2 = sigma*sigma;
sigma2_2 = 2*sigma2;
dist = 1;
T0 = .5;
T1 = -.5;     % two-sided threshold
% $$$ Ps_0 = 4*qf(dist/(2*sigma))  - 4*qf(dist/(2*sigma))^2;   % true error
% $$$ fprintf('True Ps_0=%g\n',Ps_0);
nonmixisnoise = 0;

zetalist = .5;
zetalist = 0.2:.1:1;
% zetalist = 0.9:.1:1;
dim = 4;     % number of dimensions of data
x = zeros(dim,1);  % transmitted symbol
Ps_0 = 1 - (1 - 2*qf(dist/(2*sigma)))^dim;
fprintf('True Ps_0=%g\n',Ps_0);
K = sigma2*eye(dim);

% sigma = sqrt(0.05)
% n=1 
% $$$ gammabest=0.179284  7.85382  6.52808  
% $$$ zetabest=0.2  0.5  0.5  
%n = 2
% $$$ gammabest=0.207262  7.76397  9.14323  
% $$$ zetabest=0.2  0.5  0.5  
%n = 3
% $$$ gammabest=0.180843  7.80305  11.6169  
% $$$ zetabest=0.2  0.5  0.6  
% n = 4
% $$$ gammabest=0.152615  7.81381  13.7253  
% $$$ zetabest=0.2  0.5  0.6


% sigma = sqrt(0.02);
% n = 1
% $$$ gammabest=325.785  308.182  172.636  
% $$$ zetabest=0.7  0.5  0.5 
% n = 2
% $$$ gammabest=139.186  308.834  248.086  
% $$$ zetabest=0.8  0.5  0.5  

% sigma = sqrt(0.01)
%n = 1
% $$$ gammabest=614467  154342  
% $$$ zetabest=0.5  0.5  0.5 
%n = 2
% $$$ gammabest=383446  151758  88620.8  
% $$$ zetabest=0.6  0.5  0.5 
% n = 3
% $$$ gammabest=118.159  310.284  297.543  
% $$$ zetabest=0.8  0.5  0.5  
% n = 4
% $$$ gammabest=139.574  308.415  335.304  
% $$$ zetabest=0.8  0.5  0.5  

n = 1;
nexp = 10;   % number to loop over to compute variance
N = 500; 100000;  % number of samples to estimate error


d21 = dim/2-1;

gammabest = zeros(1,3);
for zeta = zetalist
   fprintf('zeta=%g\n',zeta);
   r = zeta*zeta/sigma2_2;    % used for spherical noise
   nsqrt2r = n*sqrt(2*r);
   const1 = 1/(sqrt(2*pi)^dim * sqrt(det(K)));
   const2 = 1/(sqrt(2*pi)^dim * sqrt(det(K)))*sqrt(n)^dim * exp(-n*r)*...
         sqrt(2)^(dim/2 - 1)/sqrt(r)^(dim/2 - 1)*gamma(dim/2)*sqrt(n)^d;
   
   savePbMC = zeros(1,nexp);
   savePbIS = zeros(1,nexp);
   for nexp1 = 1:nexp
      ct = 0;    % count errors for true noise distribution
      ctis = zeros(1,3);  % count errors for IS distribution
      gammadenom = zeros(1,3);

      for k=1:N
         noiseMC = sigma*randn(dim,1);
         y = x + noiseMC;  % generate according to true noise distribution
         
         if(any(y > T0) || any(y < T1))
            ct = ct + 1;
         end

         % do different kinds of IS
         % 1: just adjust one direction
         % 2: adjust in 2 directions
         % 3: pick point on sphere

         % 1: 1 direction 
         %  noise ~ N(zeta,sigma2);
         cv = zeros(dim,1);  cv(1) = zeta;
         noiseIS = sigma*randn(dim,1) + cv;
         yis = x + noiseIS;
         if(any(yis > T0) || any(yis < T1)) % count error for sampling distr.
            % compute importance sampling weight
            w1 = exp(-1/sigma2_2*(norm(yis - x)^2));
            w2 = exp(-1/sigma2_2*(norm(yis - x - cv)^2));
            w = w1/w2;
            % fprintf('w1=%g  w2=%g  w=%g\n',w1,w2,w);
            ctis(1) = ctis(1) + w*1;
            gammadenom(1) = gammadenom(1) + w*w;
         end
         
         % 2: mixture
         % do noise that is a mixture -- bias both ways at
         % random in a random coordinate direction
         p1 = 2*(rand>0.5) - 1;      % +/- 1
         p2 = floor(dim*rand) + 1;   % pick a dimension
         mis = zeros(dim,1);
         mis(p2) = zeta*p1;
         noiseIS2 = sigma*randn(dim,1) + mis;
         yis2 = x + noiseIS2;

         if(any(yis2 > T0) || any(yis2 < T1)) % count error for sampling distr.
                                            % compute importance sampling weight
            w1 = exp(-1/sigma2_2*(norm(yis2 - x)^2));
            w2 = 0;
            cv = zeros(dim,1);
            for i1 = 1:dim
               cv(i1) = zeta;
               w2 = w2 + exp(-norm(yis2 - x - cv)^2/sigma2_2);
               cv(i1) = -zeta;
               w2 = w2 + exp(-norm(yis2 - x - cv)^2/sigma2_2);
               cv(i1) = 0;
            end
            w2 = w2/(2*dim);
            w = w1/w2;
            ctis(2) = ctis(2) + w*1;
            gammadenom(2) = gammadenom(2) + w*w;
         end

         % 3: Spherical
         u = randn(dim,1);  u = u/norm(u);   % random on unit circle
         u = zeta*u;                        % circle of radius zeta
         noiseIS3 = sigma/sqrt(n)*randn(dim,1) + u;
         yis3 = x + noiseIS3;
         if(any(yis3 > T0) || any(yis3 < T1)) % count error for sampling distr.
            % compute importance sampling weight
            yn = norm(yis3 - x);
            yn2 = yn*yn;
            ynsigma = yn/sigma;
            w1 = const1*exp(-1/sigma2_2*yn2);

            w2 = const2*exp(-1/sigma2_2*yn2*n)*...
                 besseli(d21, nsqrt2r*ynsigma)/( n*(ynsigma)^d21);
            
            w = w1/w2;
            ctis(3) = ctis(3) + w*1;
            gammadenom(3) = gammadenom(3) + w*w;
         end
         
      end  % for k
      PbMC = ct/N;
      PbIS(1) = ctis(1)/N;
      PbIS(2) = ctis(2)/N;
      PbIS(3) = ctis(3)/N;

      savePbMC(nexp1) = PbMC;
      savePbIS(1,nexp1) = PbIS(1);
      savePbIS(2,nexp1) = PbIS(2);
      savePbIS(3,nexp1) = PbIS(3);

      fprintf('   nexp1=%d   PbMC=%g  PbIS=(%g,%g,%g)\n',nexp1,PbMC,...
         PbIS(1),PbIS(2),PbIS(3));
   end % for nexp1
   varMC = var(savePbMC);
   varIS(1) = var(savePbIS(1,:));
   varIS(2) = var(savePbIS(2,:));
   varIS(3) = var(savePbIS(3,:));

   gammarat(1)=(PbIS(1) - PbIS(1)*PbIS(1))/ (gammadenom(1)/N - PbIS(1)*PbIS(1));
   gammarat(2)=(PbIS(2) - PbIS(2)*PbIS(2))/ (gammadenom(2)/N - PbIS(2)*PbIS(2));
   gammarat(3)=(PbIS(3) - PbIS(3)*PbIS(3))/ (gammadenom(3)/N - PbIS(3)*PbIS(3));
   if(gammarat(1) > gammabest(1))
      gammabest(1) = gammarat(1);
      zetabest(1) = zeta;
   end
   if(gammarat(2) > gammabest(2))
      gammabest(2) = gammarat(2);
      zetabest(2) = zeta;
   end
   if(gammarat(3) > gammabest(3))
      gammabest(3) = gammarat(3);
      zetabest(3) = zeta;
   end
   
   fprintf('   varMC=%g  varIS=',varMC);
   fprintf('%g ',varIS);
   fprintf('\n   varMC/varIS=');
   fprintf('%g ',varMC ./ varIS);
   fprintf('  gammarat=');
   fprintf('%g ',gammarat);
   fprintf('\n');

end % for zeta

fprintf('gammabest=');  fprintf('%g  ',gammabest);
fprintf('\nzetabest=');  fprintf('%g  ',zetabest);
fprintf('\n');

   


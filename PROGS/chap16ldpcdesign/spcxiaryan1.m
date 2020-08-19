
% test single-parity check decoding
%  Uses the method of Xia Ryan.  Assumes s = (1,1,1,...,1).

N = 20;
K = N-1;
M = N-K;
R = K/N;

EbN0dBlist = 20;  % 4:.5:12;   % 2:.5:4;  % 6:.5:8;
Maxiter = 2000;
s = ones(N,1);
b = ones(N,1);  b(1) = 0;  b(2) = 0;
explist = zeros(1,Maxiter);

for EbN0dB = EbN0dBlist
   fprintf('EbN0=%g   N=%d\n',EbN0dB,N);
   EbN0 = 10^(EbN0dB/10);
   sigma2 = 1/(2*R*EbN0);    % Ec = 1;
   sigma = sqrt(sigma2);
   sigma2_2 = 2*sigma2;
   nbiterr = 0;
   nbiterrIS1 = 0;
   v = 0;
   
   L = 2/sigma2;
   for niter=1:Maxiter
      noise = sigma*randn(N,1);  % noise (MC sampling)
      y = s + noise;        % received signal
      l = L*y;
      D = tanhNeq(l(1),l(2:end));
      if(D < 0)   % error event
         nbiterr = nbiterr + 1;
      end

      yIS1 = b + noise;
      lIS1 = L*yIS1;    % log likelihood ratio
      DIS1 = tanhNeq(lIS1(1),lIS1(2:end));
      if(DIS1 < 0 && all(lIS1(2) < lIS1(3:end)))
         %      if(DIS1 < 0 && lIS1(2) < lIS1(3))   % an E2 error event
         w1 = exp((yIS1(1) + yIS1(2) - 1)/sigma2);
         nbiterrIS1 = nbiterrIS1 + w1;
         v = v + w1*w1;
      end
   end

   PbMC = nbiterr / Maxiter;
   PbIS1 = nbiterrIS1 / Maxiter * (N-1);
   Pbbound = (N-1)*qf(sqrt(2)/sigma);
   % v = v/Maxiter - Pbbound*Pbbound;
   v = v/Maxiter - PbIS1*PbIS1;
   epsi = 0.1;
   % L = ceil(v/(epsi*epsi*Pbbound*Pbbound));
   L = ceil(v/(epsi*epsi*PbIS1*PbIS1));
   
   fprintf('EbN0dB=%g  sigma=%g  PbMC=%g PbIS1=%g  Pbbd=%g  L=%g\n',...
           EbN0dB,sigma,PbMC,PbIS1,Pbbound,L);
end  % for EbN0


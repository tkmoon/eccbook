
% plot the bounds on performance for a convolutianal code.
% use T(D,N) = D^5 N/(1-2DN),  so T(D) = D^5/(1-2D)
% dTdN = D^5/(1-2DN)^2   so dTdN_N=1  = D^5/(1-2D)^2
% for a R=1/2 code with dfree = 5;

% describe the code
k = 1;
n = 2;
R = k/n;
dfree = 5;
anfree = 1;  % number of paths of this length
Nfree = 1;  % number of nonzero input bits on shortest path

% plot range
EbN0dB = 4:16;

% Initialize the plot arrays
uncoded = [];							% uncoded performance
Pehard = [];							% node error for hard decoding
Pbhard = [];							% bit error for hard decoding
Pbhardlower = [];						% lower bound for hard decoding
Pesoft1 = [];							% node error for soft decoding, bound 1
Pesoft2 = [];							% node error for soft decoding, bound 2
Pbsoft1 = [];							% bit error for soft decoding, bound 1
Pbsoft2 = [];							% bit error for soft decoding,bound 2
Pbsoftlower = [];						% lower bound for soft decoding
  
for e = EbN0dB;							% loop over each SNR
  EbN0 = 10^(e/10);  % convert back from dB
  uncoded = [uncoded qf(sqrt(2*EbN0))];
  EcN0 = R*EbN0;
  pc = qf(sqrt(2*EcN0));  % determine equivalent BSC crossover probability
  % Hard decision upper bounds
  Z = sqrt(4*pc*(1-pc));
  Pehard = [Pehard (Z^5/(1-2*Z))];
  Pbhard = [Pbhard ((Z^5)/(1-2*Z)^2)/k];
  % compute Pdfree
  if(mod(dfree,2) == 0) 				% even dfree
	Pdfree = 0.5*nchoosek(dfree,dfree/2)*pc^(dfree/2)*(1-pc)^(dfree/2);
	for k1=dfree/2+1:dfree
	  Pdfree = Pdfree + nchoosek(dfree,k1)*pc^k1*(1-pc)^(dfree-k);
	end
  else									% odd dfree
	Pdfree = 0;
	for k1=(dfree+1)/2:dfree
	  Pdfree = Pdfree + nchoosek(dfree,k1)*pc^k1*(1-pc)^(dfree-k);
	end
  end
  % Hard decision lower bound
  Pbhardlower = [Pbhardlower Nfree*anfree*Pdfree/k];

  % soft decision upper bounds
  Z = exp(-EcN0);
  Pesoft1 = [Pesoft1 0.5*(Z^5/(1-2*Z))];
  Pesoft2 = [Pesoft2 exp(dfree*EcN0)*qf(sqrt(2*dfree*EcN0))*(Z^5)/(1-2*Z)];
  Pbsoft1 = [Pbsoft1 (0.5*Z^5/(1-2*Z)^2)/k];
  Pbsoft2 = [Pbsoft2 (exp(dfree*EcN0)*qf(sqrt(2*dfree*EcN0))*Z^5/(1-2*Z)^2)/k];
  % Soft decision lower bounds
  Pbsoftlower = [Pbsoftlower Nfree*anfree*qf(sqrt(2*dfree*EcN0))/k];

end
clf
semilogy(EbN0dB,uncoded);
hold on;
semilogy(EbN0dB,Pbhard,'--');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of bit error (bound)');
semilogy(EbN0dB,Pbhardlower,'--');

semilogy(EbN0dB,Pbsoft2,':');
semilogy(EbN0dB,Pbsoft1,'r:');
semilogy(EbN0dB,Pbsoftlower,':');

legend('Uncoded','Hard decision bounds', 'Soft decision bounds');

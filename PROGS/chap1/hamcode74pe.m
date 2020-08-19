% Plot the probability of decoding
% error for a Hamming code as a function of SNR

% SNR = 0:1:12;  % Eb/N0, in dB
SNR = 4:12;  % Eb/N0, in dB
EbN0 = 10.^(SNR/10);


% Stuff for Hamming code
% for (7,) code:
n = 7;
k = 4;
R = k/n;
dmin = 3;
EcN0 = EbN0*k/n;
Pnocode = qf(sqrt(2*EbN0));  % prob. of bit error with no coding
Pen = qf(sqrt(2*EcN0));      % prob of bit error for coded data, but
                             % without decoding
                             % 
Pasymptsoft = qf(sqrt(2*R*dmin*EbN0));

t = (dmin-1)/2;

% compute the weight enumerator polynomial
syms z Az;
Az = 1/(n+1)*((1+z)^n + n*(1-z)*(1-z^2)^((n-1)/2));
A = sym2poly(Az)

j1 = 0;
for p = Pen
  j1 = j1+1;

  % compute the P_k^j functions 
  for j=dmin:n
	for k1=0:t
	  P(k1+1,j) = 0;
	  for r=0:k1
		P(k1+1,j) = P(k1+1,j) + nchoosektest(j,k1-r)*nchoosektest(n-j,r)* ...
			p^(j-k1+2*r)*(1-p)^(n+j+k1-2*r);
	  end
	end
  end
  % Compute P(E), the probability of decoder error:
  PE(j1) = 0;
  for j=dmin:n
	s = 0;
	for k1=0:t
	  s = s + P(k1+1,j);
	end
	PE(j1) = PE(j1) + A(j+1)*s;
  end
end



% Compute B_j
Par = [1 1 0 1
	   1 0 1 1
	   0 1 1 1];
H = [eye(n-k) Par];
G = [Par' eye(k)];
B = zeros(n,1);
for i=1:2^k-1   % walk through all nonzero message vectors
  m = dec2bin(i,k)-'0';  % get message vector
  wm = sum(m);
  c = mod(m*G,2);
  wc = sum(c);
  B(wc) = B(wc) + wm;
  %if(wm ==3)
	%m
	%c
  %end
end

j1 = 0;
for p = Pen
  j1 = j1+1;

  
  % compute the P_k^j functions of (10-11) of Wicker
  for j=dmin:n
	for k1=0:t
	  P(k1+1,j) = 0;
	  for r=0:k1
		P(k1+1,j) = P(k1+1,j) + nchoosektest(j,k1-r)*nchoosektest(n-j,r)* ...
			p^(j-k1+2*r)*(1-p)^(n+j+k1-2*r);
	  end
	end
  end

% Compute Pb, the probability of decoder error:
  Pb(j1) = 0;
  for j=dmin:n
	s = 0;
	for k1=0:t
	  s = s + P(k1+1,j);
	end
	Pb(j1) = Pb(j1) + B(j)*s;
  end
  Pb(j1) = Pb(j1)/k;
end

clf
semilogy(SNR,Pnocode);
hold on
semilogy(SNR,Pen,'r--');
%semilogy(SNR,PE,'r:');
%semilogy(SNR,PE/k,'g:');
semilogy(SNR,Pb,'-.');

semilogy(SNR,Pasymptsoft,':');
lh = legend('BER with no coding','BER with coding','Decoded BER',...
	'Asymptotic soft-input decoding');
set(lh,'fontsize',15);
xh = xlabel('E_b/N_0 (dB)');
yh = ylabel('Probability of bit error P_b');
set(xh,'fontsize',15);
set(yh,'fontsize',15);

snrnocode = interp1(log(Pnocode),SNR,log(10^(-6)));
snrhammcode = interp1(log(Pb),SNR,log(10^(-6)));

s = 2;
set(gca,'fontsize',15);
line([snrnocode snrnocode],[10^(-6)/s 10^(-6)*s]);
line([snrhammcode snrhammcode],[10^(-6)/s 10^(-6)*s]);
line([snrhammcode snrnocode],[10^(-6) 10^(-6)]);
th = text(snrhammcode-1.9,.8*10^(-6),'Coding Gain');
set(th,'fontsize',15);
input('press return to save picture');
print -dps hamcode74pe.ps

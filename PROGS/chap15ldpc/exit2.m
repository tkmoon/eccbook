

load ldpcsim.mat

% To convert the decoder probabilities to an exit chart:
% 1) Extract the mutual information using getinfs, which calls
%    getinf, which performs numerical integration using the histogram
% 2) Stack the information up correctly according to the iteration number
%    using buildexit
% 3) Plot the corresponding information with doexitchart.
%
% Putting all this together, for the data BtoCDataFile and CtoBDatafile,
% the command is
% doexitchart(buildexit(getinfs(CtoBDataFile),getinfs(BtoCDataFile)),'r')
%

clf
subplot(2,2,1);
doexitchart(snr08_msgpass_exitpath01);  tl = title('SNR=0.8 dB');
set(tl,'fontsize',15);
subplot(2,2,2);
doexitchart(snr04_msgpass_exitpath01);  tl = title('SNR=0.4 dB');
set(tl,'fontsize',15);
subplot(2,2,3);
doexitchart(snr16_msgpass_exitpath01);  tl = title('SNR=1.6 dB');
set(tl,'fontsize',15);
subplot(2,2,4);
doexitchart(snr18_msgpass_exitpath01);  tl = title('SNR=1.8 dB');
set(tl,'fontsize',15);

print -dps exit2.ps

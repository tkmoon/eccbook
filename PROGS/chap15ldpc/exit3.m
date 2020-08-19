

load ldpcsim.mat

% To plot the trajectory:
% 1) Extract the mutual information using getinfs, which calls
%    getinf, which performs numerical integration using the histogram
% 2) Stack the information up correctly according to the iteration number
%    using buildexit
% 3) Plot the corresponding information with dotrajectory
%
% Putting all this together, for the data BtoCDataFile and CtoBDatafile,
% the command is
% dotrajectory(buildexit(getinfs(CtoBDataFile),getinfs(BtoCDataFile)),'r')
%

clf
subplot(2,2,1);
dotrajectory(snr08_msgpass_exitpath01);  tl = title('SNR=0.8 dB');
set(tl,'fontsize',15);
subplot(2,2,2);
dotrajectory(snr04_msgpass_exitpath01);  tl = title('SNR=0.4 dB');
set(tl,'fontsize',15);
subplot(2,2,3);
dotrajectory(snr16_msgpass_exitpath01);  tl = title('SNR=1.6 dB');
set(tl,'fontsize',15);
subplot(2,2,4);
dotrajectory(snr18_msgpass_exitpath01);  tl = title('SNR=1.8 dB');
set(tl,'fontsize',15);

input('press return');
print -dps exit3.ps

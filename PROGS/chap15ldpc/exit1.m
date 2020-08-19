

load ldpcsim.mat
% load the workspace that has the results from running the decoder algorithm.
% The code is 15k x 10k, and the average column weight is 2.66, so
% there are 40,000 connections
% BtoCDataFile --- itno x 40000 --- bit to check data (decoder probabilities)
% CtoBDataFile --- itno x 40000 --- check to bit data (decoder probabilities)
% infBtoC --- the information at each iteration of BtoCDataFile.
%    infBtoC = getinfs(BtoCDatafile)
% infCtoB --- similarly
% snr04_msgpass_exitpath01 -- the information ready to pass to doexitchart,
%    obtained using buildexit(getinfs(CtoBdata   ), getinfs(BtoCData   ))
%    for an SNR of 0.4 dB
% snr06_msgpass_exitpath01 -- similary for SNR=0.6 dB
% snr08_msgpass_exitpath01 -- similary for SNR=0.8 dB
% snr10msgpass_exitpath01 -- similary for SNR=1.0 dB
% snr12msgpass_exitpath01 -- similary for SNR=1.2 dB
% snr14msgpass_exitpath01 -- similary for SNR=1.4 dB
% snr16msgpass_exitpath01 -- similary for SNR=1.6 dB
% snr18msgpass_exitpath01 -- similary for SNR=1.8 dB

% plot histograms
ittodo = [1 2 11 15];
i1 = 0;
for i = ittodo
  i1 = i1+1;
  loghist(BtoCDataFile(i,:),i,i1);
end

print -dps exit1.ps
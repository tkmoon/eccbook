function [I] = getinf(histN,histX)
% where histN and histX are the output of
%  a [N,X] = hist(...) function call
% That is, N is the frequencies, X are the bin centers
% Center bin must be at Zero
% I is the mutual information


lnX = length(histX);

for k=1:lnX
    if(histN(k))
        Ia(k) = (histN(k)/sum(histN))*log2((histN(k))/(0.5*(histN(k)+histN(lnX+1-k))));
    else
        Ia(k)=0;
    end
end

for(k=1:lnX)
    Ib(k) = Ia(k) + Ia(lnX+1-k);
end

I = 0.5*sum(Ib);

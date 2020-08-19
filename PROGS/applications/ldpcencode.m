function parity = ldpcencode(info,ldpc)
% function parity = ldpcencode(info,ldpc)
%  Compute the parity bits for the information bits in info,
%  using the LDPC code structure in LDPC.
%  This works for the irregular repeat accumulate codes
%  such as for the DVB-T2 standard (and not for other codes!)
% 
%  This function returns the parity bits, from which the overall
%  systematic codeword can be formed.
%
%  The information for these IRA codes for DVB-T2 is in
%  getiradat.m, which picks the right table for the IRA code.
%  This is parsed into the correct format by ldpccodesetup1.m, 
%    which turns the data into a parity check matrix H
%    and write it to a file.
%
%  The file is read by readldpcfile.m.  
%    This function assumes that the index offset is 1 
%    (as set when the file is read).
%    readldpcfile produces the sparse representation of the 
%    H matrix used by the decoder.


% TKM 6/23/17

Nminusk = ldpc.M;

parity = zeros(1,Nminusk);
for m = 1:Nminusk  % for each parity
   if(m == 1)
      parity(m) = mod(sum(info(ldpc.Nm{m}(1:end-1))),2);
   else
      parity(m) = mod(sum(info(ldpc.Nm{m}(1:end-2))),2);
   end
end
parity = mod(cumsum(parity),2);


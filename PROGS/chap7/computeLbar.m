function L = computeLbar(n,k,q,t)
%function L = computeLbar(n,k,q,t)
%
% Compute the average number of codewords in a
% Hamming sphere of radius t around an randomly chosen
% point for an (n,k) code over GF(q)

% Todd K. Moon, Feb. 12, 2004

L = 0;
for s=0:t
  L = L + nchoosek(n,s)*(q-1)^s;
end
L = L/(q^(n-k));


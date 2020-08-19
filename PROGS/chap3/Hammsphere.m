function V = Hammsphere(n,q,t)
%function V = Hammsphere(n,q,t)
%
% Compute the number of vectors in a Hamming sphere of radius t
% in a q-ary code of length n

% Todd K. Moon, May 25, 2004

V = 0;
for s=0:t
  V = V + nchoosek(n,s)*(q-1)^s;
end


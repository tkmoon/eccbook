function nk = nchoosektest(n,k)


if(k < 0 | k > n)
  nk = 0;
else
  nk = nchoosek(n,k);
end

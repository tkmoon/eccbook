int computetm(int n,int k,int m)
// Compute tm, the number of symbols that can be corrected, for a 
// (n,k) Reed-Solomon with the Guruswami(m) decoding algorithm
//
{
   int Km, tm;
   int C;
   int v = k-1;
   int lk;
   int n1;

   if(m==0) {
	  Km = int(ceil((n+v+1.)/2.));
	  tm = n-Km;
	  return tm;
   }

   C = n*m*(m+1)/2;
   int l,k1 = 0;
   while(1) {
	  k1 = k1+1;
	  l = (m*k1)-1;
	  lk = int(floor(double(l)/v));
	  n1 = (l+1)*(lk+1) - v*lk*(lk+1)/2;
	  if(n1 > C) break;
   }
   Km = k1;  
   tm = n-Km;
   return tm;

}

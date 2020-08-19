int computeLm(int n,int k, int m)
// Lm = computeLm(n,k,m);
//
// Compute the maximum length of a list for an (n,k)
// code using GS(m) decoding

// Todd K. Moon, Feb. 12, 2004
{
   double v = k-1;
   double t;
   if(m==0) return 1;
   t = (v+2)/(2*v);
   int Lm = int(floor(sqrt( n*m*(m+1)/v + t*t) - t));
   return int(Lm);
}

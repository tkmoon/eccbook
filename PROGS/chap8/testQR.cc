//  Program: testQR.cc --- do the arithmetic for a QR code
//  Todd K. Moon

#include "GFNUM2m.h"
#include "polynomialT.cc"

// create instantiations of the polynomial class of type GFNUM2m
template class polynomialT<GFNUM2m>;


int main()
{
   int i;
   int p = 23;
   int s = 2;
   // determine the field size
   int m;
   for(m = 2; m < 20; m++) {
	  if(((1<<m)-1) % p == 0) {
		 cout << "m=" << m << endl;
		 break;
	  }
   }
   GFNUM2m::initgf(m);			// initialize field
   // find the quadratic residues
   int *Q = new int[(p-1)/2];
   int j;
   cout << "Quadratic residues: ";
   for(i = 1,j=0; i <= p/2; i++,j++) {
	  Q[j] = (i*i) % p;
	  cout << Q[j] << " ";
   }
   cout << endl;

   POLYC(GFNUM2m,l1,{0,1});		// linear term x + 0
   polynomialT<GFNUM2m> q(1);	// generator polynomial q(x)
   polynomialT<GFNUM2m> n(1);	// generator polynomial n(x)

   int inset(int, int *, int);

   GFNUM2m beta = A^(((1<<m)-1)/p);
   cout << "beta=" << beta << endl;
   for(i = 1; i <p; i++) {
	  l1[0] = beta^i;
	  if(inset(i,Q,(p-1)/2)) {		// i is a quadratic residue
		 q *= l1;
	  }
	  else {					// i is not a quadratic residue
		 n *= l1;
	  }
   }
   cout << "q=" << q << endl;
   cout << "n=" << n << endl;
}


int inset(int m, int *S, int n)
// see if m is in S, where S has n elements
{
   int i;
   for(i = 0; i < n; i++) {
	  if(S[i] == m) return 1;
   }
   return 0;
}
   

/*
Local Variables:
compile-command: "g++ -o testQR -g testQR.cc GFNUM2m.cc"
End:
*/



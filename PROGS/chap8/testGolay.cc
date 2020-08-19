//  Program: testGolay.cc --- test the algebraic Golay decoder
//  Todd K. Moon, March 30, 2004 

#include "GFNUM2m.h"
#include "polynomialT.cc"

// create instantiations of the polynomial class of type GFNUM2m
template class polynomialT<GFNUM2m>;
template class polytemp<GFNUM2m>;

int main()
{
   int i,j,l;
   int j1,j2,j3;
   GFNUM2m::initgf(11);			// initialize the class
   int n = 23;
   polynomialT<GFNUM2m> r;
   r.buildspace(22);			// make space for the data
   polynomialT<GFNUM2m> L;		// error locator polynomial
   L.buildspace(3);				// make space for the data
   L[3] = 1;					// set high-order coefficient of poly.

   GFNUM2m s1,s3,s9;			// syndromes
   GFNUM2m sigma1, sigma2, sigma3; // coefficients of error locator
   GFNUM2m beta;				// 23rd primitive root of unity
   GFNUM2m D;					// cubic term
   GFNUM2m t;					// temporary
   GFNUM2m x;					// cube root of D

   beta = A^89;
   int nroots, degree;

   for(j = 0; j < n; j++) r[j] = 0;  // clear out previous

   // walk over all patterns of 1, 2 or 3 errors
   for(j1 = 0; j1 < n; j1++) {
	  cout << "j1=" << j1 << endl;
	  for(j2 = 0; j2 < n; j2++) {
		 for(j3 = 0; j3 < n; j3++) {
			r[j1] = r[j2] = r[j3] = 1;

			// rather than call a function, we'll just do all the 
			// work right here
			
			// evaluate the three syndromes first
			s1 = r(beta);
			s3 = r(beta^3);
			s9 = r(beta^9);
			if(s1 == 0 && s3 == 0 && s9 == 0) {
			   cout << "no errors\n";
			}
			else if((s1^9) == s9 && (s3^3) == s9) { // single error 
			   L[3] = 0;
			   L[2] = 0;
			   L[1] = 1;
			   L[0] = s1;
			   degree = 1;
			}
			else {				// two or three errors
			   t = (s1^3) + s3;
			   D = t*t + ((s1^9) + s9)/t;
			   sigma1 = s1;
			   x = D^1365;
			   sigma2 = s1*s1 + x;
			   sigma3 = s3 + s1*x;
			   L[0] = sigma3;
			   L[1] = sigma2;
			   L[2] = sigma1;
			   L[3] = 1;
			   degree = 3;
			}

			// now find the roots.  Rather than use a chien searh,
			// we'll just sweep through the coefficients
			int nroots = 0;
			for(i = 0; i < n; i++) {
			   if(L(beta^i) == 0) { // root found
				  r[i] = 0;
				  nroots++;
				  if(nroots == degree) { // done!
					 break;
				  }
			   }
			}
			// Decoding done!  

			// Now check for remaining errors
			for(i = 0; i < n; i++) {
			   if(r[i] != 0) {
				  cout << "Error: j1=" << j1 << " j2=" << j2 << 
					 " j3=" << j3 << "  i=" << i << endl;
			   }
			}
			r[j1] = r[j2] = r[j3] = 0;  // clear out to get ready again
		 }			   
	  }
   }
}

/*
Local Variables:
compile-command: "g++ -o testGolay -g testGolay.cc GFNUM2m.cc"
End:
*/



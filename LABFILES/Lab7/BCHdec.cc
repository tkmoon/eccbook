// BCHdec.cc -- a general BCH decoder
// Todd K. Moon
#include "BCHdec.h"
//#include "berlmass2.cc"
#include "berlmassBCH2.cc"

//#include "berlmassBCH.cc"
//#include "polynomialT.cc"

//#include "gcdpoly.cc"

BCHdec::BCHdec(int t_in, int n_in, GFNUM2m A1_in) :
   Searcher(t_in)
{
   t = t_in;
   t2 = 2*t_in;
   n = n_in;
   A1 = A1_in;
   s = new GFNUM2m[2*t];
   Lambda = new GFNUM2m[t+1];
}

void
BCHdec::decode(GF2 *r, GF2 *dec)
{
   int i,j;
   int errloc;
   GFNUM2m p;
   GFNUM2m x;
   GFNUM2m *rs;					// array of roots

   // Step 1: evaluate the syndromes

   // Fill in the Blanks...

   // Step 2: Determine the error locator polynomial
   berlmassBCH2(s, t2,  Lambda, nu);
   // (or berlmass or berlmass2...)

   // Step 3: Find the roots of the error locator polynomial
   int nroots; 
   rs = Searcher.Search(Lambda,nu,nroots); // from ChienSearch

   // Step 4: Find error values: Not necessary for binary BCH codes

   // Step 5: Correct errors corresponding to inverse roots
   for(i = 0; i < n; i++) {  // copy over the input bits
	  dec[i] = r[i];

   }
   for(i = 0; i < nroots; i++) {
	  errloc = (GFNUM2m::getN() - rs[i].getp()) % GFNUM2m::getN();
	  dec[errloc] = 1+dec[errloc];  // complement the bits
   }
}


/*
Local Variables:
compile-command: "g++ -c -g BCHdec.cc"
End:
*/
	  

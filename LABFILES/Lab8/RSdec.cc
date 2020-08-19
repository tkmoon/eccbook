// RSdec.cc -- a general Reed-Solomon decoder
// Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "RSdec.h"
#include "berlmass2.cc"

#include "berlmass.cc"

//#include "gcdpoly.cc"
// #define BMALG					// define this to use BM alg. for decoding
// otherwise, Euclidean algorithm is used for decoding


// construtor
RSdec::RSdec(int t_in, int n_in, GFNUM2m A1_in)
   :Searcher(t_in)
{
   t = t_in;
   t2 = 2*t_in;
   n = n_in;
   A1 = A1_in;
   s = new GFNUM2m[2*t];
   Lambda = new GFNUM2m[t+1];
   spoly.setc(t2,1);  // allocate space for this
   Lambdap.setc(t,1); // allocate space for this
}

int RSdec::decode(unsigned char *r, unsigned char *dec)
{
   int i;
   GFNUM2m rgf[n];
   GFNUM2m decgf[n];
   for(i = 0; i < n; i++) {
	  rgf[i] = r[i];
   }
   int ret = decode(rgf,decgf);
   for(i = 0; i < n; i++) {
	  dec[i] = decgf[i].getv();
   }
   return ret;
}


#ifdef BMALG
int
RSdec::decode(GFNUM2m *r, GFNUM2m *dec)
{
   int i,j;
   int errloc;
   GFNUM2m p,dp2;			  // fast computation of formal derivative
   GFNUM2m x;
   GFNUM2m *rs;					// array of roots

   for(i = 0; i < n; i++) {  // copy over the input bits
	  dec[i] = r[i];
   }
   // Step 1: evaluate the syndromes
   // fill in the blanks ...

   // Step 2: Determine the error locator polynomial
   // fill in the blanks ...

   // Step 3: Find the roots of the error locator polynomial
   // fill in the blanks ...

   // Step 4: Find error values: Not necessary for binary BCH codes
   // fill in the blanks ...
   return 1;
}


#else
#include "gcdpoly.cc"

// decoding using the Euclidean algorithm
int
RSdec::decode(GFNUM2m *r, GFNUM2m *dec)
{
   int i,j;
   int errloc;
   GFNUM2m p,dp;
   GFNUM2m dp2;					// fast computation of formal derivative
   GFNUM2m x;
   GFNUM2m *rs;					// array of roots

   for(i = 0; i < n; i++) {  // copy over the input bits
	  dec[i] = r[i];

   }
   polytemp<GFNUM2m> *sS = polytemp<GFNUM2m>::gettemppoly(2*t-1);
   // Step 1: evaluate the syndromes
   // fill in the blanks ...

   // Step 2: Solve the key equation
   // fill in the blanks ...

   // Step 3: Find the roots of the error locator polynomial
   return 1;
}


#endif
/*
Local Variables:
compile-command: "g++ -c -g RSdec.cc"
End:
*/
	  

// GFNUM2m.cc --- Function definitions for GFNU2m,
// performing Galois field arithmetic
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "GFNUM2m.h"

// Define static variables in GFNUM2m
int *GFNUM2m::p2v = 0;			// convert exponent to vector
int *GFNUM2m::v2p = 0;			// a list of elements to convert from
								// vector to exponential notation
int GFNUM2m::gfm = 0;			// vector size of field element
int GFNUM2m::gfN = 0;			// number of nonzero-elements in field
outformat GFNUM2m::outtype = GFpower;
								// default to exponential output
// ALPHA elements from the field
GFNUM2m ALPHA;                  // define the element alpha
GFNUM2m& A= ALPHA;              // and a shorthand reference to it


void GFNUM2m::initgf(int m)
{
   //do the initialization by using only the size of the field
   // m in GF(2^m).
   // A fixed set of primitive polynomials is used.
   // Octal:
   unsigned int g[] = {1,1,7,013, 023, 045, 0103, 0211, 0435, 01021, 02011,
					   04005, 010123, 020033, 042103, 0100003};
   if(m>sizeof(g)/sizeof(unsigned int)) {
	  cerr << "Error: must specify connection polynomial for m" << endl;
   }
   GFNUM2m::initgf(m,g[m]);
}

// initgf: (1) Build the v2p and p2v tables
//         (2) set the global variable ALPHA
//         (3) set the static member variables gfm and gfN
void GFNUM2m::initgf(int m,unsigned int g)
// m -- GF(2^m)
// g -- generator polynomial, bits represent the coefficients:
// e.g. g = 0x13 = 1 0011 = D^4 + D + 1
{
   int i,j;


   if(m > sizeof(unsigned int)*8) {	// too many bits!
	  cerr << "Error: Degree too large in GFNUM2m" << endl;
   }
	  
   ALPHA.v = 2;					// set up alpha element
   gfm = m;
   gfN = (1<<m)-1;				// gfN = number of nonzero field elements,
								// gfN = 2^n -1

   if(v2p) delete[] v2p;		// delete any prior stuff
   if(p2v) delete[] p2v;
   
   v2p = new int[gfN+1];  		// table to convert vector to power form
   p2v = new int[gfN+1];		// tabel to convert power to vector form

   // convert g to an integer
   int gint=0;
   int mask1 = 1<<(m-1);
   int bitout;
   gint = g;
   v2p[0] = gfN;				// no conversion of this element

   int lfsr = 1;				// initial load of lfsr representation: alpha^0=1
   for(i = 1; i <= gfN; i++) {
	  v2p[lfsr] = i-1;
	  bitout = ((lfsr&mask1)!=0);
	  lfsr = ((lfsr<<1) ^ bitout*gint)&gfN;
   }

   p2v[0] = 1;
   for(i = 1;i <= gfN; i++) {
	  p2v[v2p[i]] = i;
   }
   p2v[gfN] = 1;
}


ostream& operator<<(ostream& s,const GFNUM2m& arg)
{
   if(arg.getouttype() == GFpower) {
	  if(arg.v < 2) return s << arg.v;
	  else {
		 int e = arg.v2p[arg.v];
		 if(e == 1) return s << "A";
		 else return s << "A^" << e;
	  }
   }
   else {						// vector (numeric) output
	  return s << arg.v;
   }
}


/*
Local Variables:
compile-command: "g++ -c -g GFNUM2m.cc"
End:
*/

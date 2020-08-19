// RSenc.h -- a general RS encoder
// Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef RSenc_H
#define RSenc_H

#include "GFNUM2m.h"
#include "polynomialT.h"

class RSenc {
   int n;						// block length
   int k;						// number of input symbols
   int t; 						// random error correcting capability
   int j0;						// exponent (j0=1 for narrow sense)
   GFNUM2m A1;					// primitive element used
   polynomialT<GFNUM2m> g;		// generator polynomial
   polynomialT<GFNUM2m> m;		// message polynomial
   polynomialT<GFNUM2m> c;		// code polynomial
public:
   RSenc(int n_in, int k_in, int t_in, int j0=1, GFNUM2m A1_in = ALPHA);
   ~RSenc() { };
   void encode(unsigned char *m_in, unsigned char *c_out); 
             // pass in array of data to be encoded, return array of encoded
   GFNUM2m *encode(GFNUM2m *m_in); 
             // pass in array of data to be encoded, return pointer to array
};
#endif

/*
Local Variables:
compile-command: "g++ -c -g RSenc.cc"
End:
*/


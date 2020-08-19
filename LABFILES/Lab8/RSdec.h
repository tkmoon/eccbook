// RSdec.h -- a general RS decoder
// Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef RSdec_H
#define RSdec_H

#include "GFNUM2m.h"
#include "polynomialT.h"
#include "ChienSearch.h"

class RSdec {
   int t; 						// error correcting capability
   int t2;						// 2*t
   int n;						// block length
   int nu;						// degree of error locator polynomial
   GFNUM2m A1;					// primitive element used
   ChienSearch Searcher;		// Chien search object
   GFNUM2m *s;					// syndrome storage
   GFNUM2m *Lambda;				// error locator polynomial
   polynomialT<GFNUM2m> spoly;	// syndrome polynomial
   polynomialT<GFNUM2m> Lambdap;// lambda polynomial
   polynomialT<GFNUM2m> Omega;	// error evaluator polynomial
public:
   RSdec(int t_in, int n_in, GFNUM2m A1_in = ALPHA);
   ~RSdec() {  delete[] s; delete[] Lambda;};
   int decode(GFNUM2m *r, GFNUM2m *dec);
   // given a received vector r, return the decoded
   // vector in dec
   int decode(unsigned char *r, unsigned char *dec);
   // given a received vector r (whose elements are 
   // characters to be interpreted as field elements for
   // fields up to GF(2^8), return the decoded
   // vector in dec
};



#endif

/*
Local Variables:
compile-command: "g++ -c -g RSdec.cc"
End:
*/


// BPSKmodvec.h --- a simple BPSK modulator to modulate vectors of data
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef BPSKMODVEC_H
#define  BPSKMODVEC_H
#include "BPSKmod.h"

class BPSKmodvec 
{
   int n;						// length of vector
public:
   double *v;					// vector of data
   BPSKmod mod1;				// scalar modulator
   BPSKmodvec(int in_n,double a_in=1) : mod1(a_in) {
	  n = in_n;
	  v = new double[n];
   }
   ~BPSKmodvec() {
	  delete[] v;
   }
   double *mod(unsigned int *bitsin) {
	  for(int i = 0; i < n; i++) {
		 v[i] = mod1.mod(bitsin[i]);
	  }
	  return v;
   }
};


#endif

/*
Local Variables:
compile-command: "g++ -c BPSKmodvec.h"
End:
*/


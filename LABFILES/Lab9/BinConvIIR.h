// BinConvIIR.h -- declarations for an (n,k) polynomial
// binary convolutional coder
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef BINCONVIIR_H
#define BINCONVIIR_H
#include "BinConv.h"
#include "BinTF2.h"
#include "matalloc.h"			// matrix allocation stuff

// The code is represented as a systematic transfer function matrix
// G = [1   0 ... 0  h1n/g
//      0   1 ... 0  h2n/g
//      ...
//      0   0 ... 1 hkn/g]
// 
// The last column of G consists of transfer functions with the
// same denominator
// The polynomial description is the same as in TF2


class BinConvIIR : public BinConv {
   BinTF2 tf;					// multi-input transfer function
public:
   BinConvIIR(int in_k, int in_n, int degs, unsigned int* h_in,
			  unsigned int g_in);
   ~BinConvIIR() {  };
   unsigned int *encode(const unsigned char *ins); // encode one step of input
   unsigned int getstate() const;		// return the state of the encoder
   void setstate(unsigned int state); // set the state of the encoder
};

#endif

/*
Local Variables:
compile-command: "gcc -c -g BinConvIIR.cc"
End:
*/

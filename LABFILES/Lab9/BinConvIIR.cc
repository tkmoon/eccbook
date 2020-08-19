// BinConv.cc -- (n,k) binary convolutional coder
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only


#include "BinConvIIR.h"			// the multi-input convolutional object
#include "matalloc.h"
#include <iostream>
using namespace std;

BinConvIIR::BinConvIIR(int in_k, int in_n, int in_p, unsigned int* h_in,
				 unsigned int g_in)
   : BinConv(in_k,in_n), tf(in_p,in_k,h_in,g_in)
{
   nu = in_p;
};


// encode ins[0] ... ins[k-1] to get the n outputs
unsigned int *
BinConvIIR::encode(const unsigned char *ins)
{
   // fill in the blanks ...

   return outs;
}

unsigned int
BinConvIIR::getstate() const
{
   // fill in the blanks ..

   return 0;  // return the actual state...
}

void
BinConvIIR::setstate(unsigned int instate)
{

   // fill in the blanks ...
}

/*
Local Variables:
compile-command: "g++ -c -g BinConvIIR.cc"
End:
*/

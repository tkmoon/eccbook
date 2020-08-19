// BinConv.cc -- (n,k) binary convolutional coder
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "BinConvFIR.h"			// the multi-input convolutional object
#include "matalloc.h"
#include <iostream>
using namespace std;

BinConvFIR::BinConvFIR(int in_k, int in_n, int *degs, unsigned int** h_in)
   : BinConv(in_k,in_n)
{
   int i,j;

   CALLOCMATRIX(h,unsigned int, k,n);
   nui = new int[k];
   mem = new unsigned int[k];
   maxdeg = 0;
   nu = 0;
   for(i = 0; i < k; i++) {
	  nui[i] = degs[i];
	  if(nui[i] > maxdeg) maxdeg = nui[i];
	  nu += nui[i];
	  for(j = 0; j < n; j++) {
		 h[i][j] = h_in[i][j];
	  }
   }
}


// encode ins[0] ... ins[k-1] to get the n outputs
unsigned int *
BinConvFIR::encode(const unsigned char *ins)
{
   // Fill in the Blanks ...
   return outs;
}

unsigned int
BinConvFIR::getstate() const
{
   // fill in the blanks ...
   
   return 0;    // return the actual state ...
}

void
BinConvFIR::setstate(unsigned int state)
{  
   // Fill in the blanks ...
}

/*
Local Variables:
compile-command: "g++ -o testconvFIR -g testconvFIR.cc BinConvFIR.cc"
End:
*/

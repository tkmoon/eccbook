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
std::cout << "in BinConvFIR constructor" << std::endl;
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
std::cout << "done with BinConvFIR constructor" << std::endl;
}


// encode ins[0] ... ins[k-1] to get the n outputs
unsigned int *
BinConvFIR::encode(const unsigned char *ins)
{
   // Fill in the blanks ...

   return outs;
}

unsigned int
BinConvFIR::getstate() const
{
   int i;
   unsigned int state = 0;
   int p = 0;
   for(i = 0; i < k; i++) {
	  state |= ((mem[i]&((1<<nui[i])-1))<<p);
	  p += nui[i];
   }
   return state & ((1<<nu)-1);
}

void
BinConvFIR::setstate(unsigned int state)
{  
   int i;
   unsigned int st;
   for(i = 0; i < k; i++) {
	  st = (state & ((1<<nui[i])-1)); // mask off piece of state
	  mem[i] = st;
	  state = state >> nui[i];
   }
}

/*
Local Variables:
compile-command: "g++ -o testconvFIR -g testconvFIR.cc BinConvFIR.cc"
End:
*/

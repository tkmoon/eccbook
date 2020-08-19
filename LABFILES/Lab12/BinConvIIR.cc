// BinConv.cc -- (n,k) binary convolutional coder
// Todd K. Moon

#include "BinConvIIR.h"			// the multi-input convolutional object
#include "matalloc.h"
#include <iostream>
using namespace std;

BinConvIIR::BinConvIIR(int in_k, int in_n, int in_p, unsigned int* h_in,
				 unsigned int g_in)
   : BinConv(in_k,in_n), tf(in_p,in_k,h_in,g_in)
{
   nu = in_p;
//   prevstate = 0;
//   outputmat = 0;
};


// encode ins[0] ... ins[k-1] to get the n outputs
unsigned int *
BinConvIIR::encode(const unsigned char *ins)
{

   // Fillin the blanks ...
   return outs;
}

unsigned int
BinConvIIR::getstate() const
{
   unsigned int state, flipstate=0;
   int i;
   state =  tf.getstate();
   // since TF puts things in the nonconventional direction, 
   // reverse the order of the bits
   for(i = 0; i < nu; i++) {
	  if(state & (1<<i)) flipstate |= (1<<(nu-i-1));
   }
   return flipstate;

//   return state;
}

void
BinConvIIR::setstate(unsigned int instate)
{
   unsigned int flipstate = 0;
   // since TF puts things in the nonconventional direction, 
   // reverse the order of the bits
   int i;
   for(i = 0; i < nu; i++) {
	  if(instate & (1<<i)) flipstate |= (1<<(nu-i-1));
   }
   tf.setstate(flipstate);

//   tf.setstate(instate);
}

/*
Local Variables:
compile-command: "g++ -c -g BinConvIIR.cc"
End:
*/

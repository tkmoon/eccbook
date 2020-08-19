// BPSKmod.h --- a simple BPSK modulator to modulate a scalar of data
// Todd K. Moon

#ifndef BPSKMOD_H
#define  BPSKMOD_H

class BPSKmod 
{
public:
   double a;					// signal amplitude
   BPSKmod(double a_in=1) {
	  a = a_in;
   }
   ~BPSKmod() {}
   double mod(unsigned int bitin) {
	  double v = a*(2.0*int(bitin)-1.0);
	  return v;
   }
};


#endif

/*
Local Variables:
compile-command: "g++ -c BPSKmod.h"
End:
*/


// Turboenc.h -- a Turbo encoder
// Todd K. Moon
#ifndef TURBOENC_H
#define TURBOENC_H

#include "BinConvIIR.h"   // the systematic convolutional encoder
#include "interleave.h"   // the random interleaver

class Turboenc {
protected:
   unsigned char *inbits;		// input bits
   unsigned char *interleavebits; // permuted bits
   unsigned int *outputbits;	// the encoded output bits
   int blocklen;				// the number of input bits in coded block
   unsigned char **P;			// puncture matrix
   int puncturelen;				// width of puncture matrix
   int puncturecycle;			// column counter for puncture matrix
public:
   Turboenc(int deg, unsigned int h_in,unsigned int g_in,
			int in_blocklen,
			unsigned int interleaveseed = 0,
			unsigned char **punctureP=0, int puncturelen=0);

   ~Turboenc() { delete[] inbits; delete[]interleavebits; delete[]outputbits;}
   unsigned int *encode(const unsigned char *ins);

   BinConvIIR Enc1;				// first encoder
   BinConvIIR Enc2;				// second encoder
   interleave interleaver;		// interleaver
   double R;					// rate of code(default=1/3 without puncturing)
   unsigned int getstate1() { return Enc1.getstate(); }
   unsigned int getstate2() { return Enc2.getstate(); }

};



#endif
/*
Local Variables:
compile-command: "g++ -c Turboenc.cc"
End:
*/


// Turbodec.h -- a Turbo decoder
// Todd K. Moon
#ifndef TURBODEC_H
#define TURBODEC_H

#include "BinConvIIR.h"   // the systematic convolutional encoder
#include "interleave.h"   // the random interleaver
#include "Turboenc.h"			// the turbo encoder
#include "BCJR.h"				// the BCJR object

class Turbodec {
   Turboenc enc;				// encoder we are working on
   BCJR bcjr;					// BCJR object
   int blocklen;				// length of data block
   double **r1;					// data for encoder 1
   double **r2;					// data for encoder 2
   double **prior1, **prior2;	// locations for prior probabilities
   int numbranch;				// number of branches = 2^k
   unsigned char **P;			// puncture matrix
   int puncturelen;				// width of puncture matrix
   int puncturecycle;			// column counter for puncture matrix
   double R;					// code rate, default=1/3 with no puncturing
   double sigma2;				// noise variance
public:
   Turbodec(int deg, unsigned int h_in,unsigned int g_in,
			int in_blocklen,
			unsigned int interleaveseed = 0,
			unsigned char **punctureP=0, int puncturelen=0);
   ~Turbodec() { 
	  FREEMATRIX(prior1); FREEMATRIX(prior2);
      FREEMATRIX(r1);  FREEMATRIX(r2);
   };
   double **decode(const double *r,int numit, 
			 unsigned int finalstate1=0, unsigned int finalstate2=-1);
   // turbo decoder function.  Returns a pointer to the array
   // of posterior probabilities

   void setsigma2(double in_sigma2) { 
	  sigma2 = in_sigma2;
	  bcjr.setsigma2(in_sigma2);
   };
};



#endif
/*
Local Variables:
compile-command: "g++ -c Turbodec.cc"
End:
*/


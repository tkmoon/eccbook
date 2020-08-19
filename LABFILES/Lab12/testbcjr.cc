//  Program: testbcjr.cc
//
//  Todd K. Moon
//

#include <iostream>
using namespace std;
#include <math.h>

#include "BinConvIIR.h"
#include "matalloc.h"
#include "BPSKmodvec.h"
#include "BCJR.h"

int main()
{
   int i,j;
   int decout;

   int interleaveseed = 1;
   int k = 1;
   int n = 2;
   unsigned int gnum = 4;		// 100 = 1
   unsigned int gden = 5;		// 101 = 1+D^2
   int p = 2;
   BinConvIIR enc(k,n,p,&gnum,gden);  // build the encoder object

   unsigned char d[] = {1,1,0,0,1,0,1,0,1,1}; // input data stream
   int N = sizeof(d)/k;			// number of time steps in test array
   int cblocklen = N*n/k;		// coded block length
   unsigned int *out;
   unsigned int *encout = new unsigned int[cblocklen];
								// array for all outputs
   cout << "State sequence: ";
   for(i = 0; i < N; i++) {
	  out = enc.encode(&d[i]);
	  cout << enc.getstate() << " ";
	  encout[2*i] = out[0];
	  encout[2*i+1] = out[1];
   }
   cout << endl;
   cout << "Binary output sequence: ";
   for(i = 0; i < N; i++) {
	  cout << int(encout[2*i]) << int(encout[2*i+1]) << " ";
   }
   unsigned int finalstate = enc.getstate();
   cout << "Final state=" << finalstate << endl;
   BPSKmodvec modulator(cblocklen);
   double *modouts;
   modouts = modulator.mod(encout);	// modulate the block of bits
   // add the noise --- in this case, fixed for consistency with example
   // but this was generated originally using sigma^2 = 0.45
   double sigma2 = 0.45;
   // use the values in the example, and subtract off the modulated
   // data to get the noise
   double noise[] = {2.53008-1,0.731636-1,-0.523916-1,1.93052-1,
					 -0.793262+1,0.307327-1,-1.24029+1, 0.784426-1,
					 1.83461-1, -0.968171+1,-0.433259+1,1.26344-1,
					 1.31717-1,0.995695-1,-1.50301+1, 2.04413-1,
					 1.60015-1, -1.15293+1,0.108878-1,-1.57889+1};
   double **chanouts;			// array for channel outputs
   CALLOCMATRIX(chanouts,double,N,2);
   for(i = 0; i < N; i++) {
	  chanouts[i][0] = modouts[2*i] + noise[2*i];
	  chanouts[i][1] = modouts[2*i+1] + noise[2*i+1];
   }
   cout << "Outputs with noise: " << endl;
   MATDUMP(chanouts,N,2);

   // Declare the BCJR object
   BCJR bcjr(enc,N,sigma2);
   double **priors;				// array of prior probabilities
   CALLOCMATRIX(priors,double,N,2);
   for(i = 0; i < N; i++) { 	// set up uniform priors
	  priors[i][0] = priors[i][1] = 0.5;
   }

   // Do the forward-backward computation
   bcjr.MAP1((const double **)chanouts,(const double **)priors,finalstate);
   MATDUMP(bcjr.alpha,N,4);
   MATDUMP(bcjr.beta,N+1,4);
}


/*
Local Variables:
compile-command: "g++ -o testbcjr -g testbcjr.cc BCJR.cc BinConvIIR.cc "
End:
*/



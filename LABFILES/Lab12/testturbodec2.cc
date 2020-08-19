//  Program: testturbodec2.cc
//
//  Todd K. Moon
//

#include <iostream>
using namespace std;
#include <math.h>

#include "BinConvIIR.h"
// #include "matalloc.h"
#include "BPSKmodvec.h"
#include "Turboenc.h"
#include "Turbodec.h"
double gran(void);
double uran(void);

int main()
{
   int i,j;
   int decout;

   int k = 1;					// number of input bits
   unsigned int gnum = 0x11;   // 1 0001
   unsigned int gden = 0x1F;     // 1 1111
   int p = 4;					// degree of denominator
   int N = 65536;				// block length

   unsigned char **P=0;			// puncture matrix

   CALLOCMATRIX(P,unsigned char,2,2);  // use this to build the puncture matrix
   P[0][0] = 1;  P[0][1] = 0;
   P[1][0] = 0;  P[1][1] = 1;

   int interleaveseed = 1;

   Turboenc encoder(p,gnum,gden,N,interleaveseed,P,2);
   cout << "Encoder built" << endl;
   int cblocklen = int(N/encoder.R);
   Turbodec decoder(p,gnum,gden,N,interleaveseed,P,2);
   cout << "Decoder built" << endl;

   cout << "cblocklen=" << cblocklen << endl;
   BPSKmodvec modulator(cblocklen);
   unsigned char *out;

   unsigned char *bitin = new unsigned char[N]; // encoder inputs
   unsigned char *bitout = new unsigned char[N]; // decoder outputs
   unsigned int *encout;		// encoder output
   double *allouts;				// modulated data
   double **post_prob;			// posterior probabilities

   double SNRstart, SNRend, SNRstep, SNRdb,SNR;
   unsigned int finalstate1, finalstate2;
   unsigned long int numbits;
   unsigned long int biterrs;
   double sigma2,sigma;
   int numerrstocount = 100;
   int numit;
   
   cout << "numerrstocount=" << numerrstocount << endl;
//   numit = 1;  double SNRlist[] = {.25};  //, 1, 1.5, 2, 2.5};
//   numit = 2;  double SNRlist[] = {.25};  // 0.5,1,1.5,2,2.5};
//   numit = 3;  double SNRlist[] = {2};
//   numit = 6;  double SNRlist[] = {.25};  // 0.5,0.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5};
//   numit = 10;  double SNRlist[] = {.25};  // 0.5,0.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5};
   numit = 18;  double SNRlist[] = {.25};  // 0.5,0.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5};
   int nsnr = sizeof(SNRlist)/sizeof(double);

   int nblock = 128;
   int blockcount;
   cout << "numit=" << numit << endl;
   for(int isnr = 0; isnr < nsnr; isnr++) {
	  SNRdb = SNRlist[isnr];
	  SNR = pow(10.,SNRdb/10);
	  numbits = 0;
	  biterrs = 0;
	  sigma2 = 1/(2*SNR*encoder.R);
	  cout << "SNR(dB)=" << SNRdb << "  SNR=" << SNR << "   sigma2=" 
		   << sigma2 << endl;
	  decoder.setsigma2(sigma2);
	  sigma = sqrt(sigma2);
	  blockcount = 0;
	  do {
		 for(i = 0; i < N; i++) {	// generate the random bits
			bitin[i] = uran()>0.5;
		 }
		 numbits += N;
		 encout = encoder.encode(bitin); // encode the data
		 finalstate1 = encoder.getstate1();
		 finalstate2 = encoder.getstate2();
		 allouts = modulator.mod(encout);
		 // add the noise
		 for(i = 0; i < cblocklen; i++) {
			allouts[i] += sigma*gran();
		 }
		 post_prob = decoder.decode(allouts,numit,finalstate1,finalstate2);
		 for(i = 0; i < N; i++) { // count up the errors
			bitout[i] = post_prob[i][1] > 0.5;
			if(bitout[i] != bitin[i]) {
			   biterrs++;
			}
		 }
		 cout << "Biterrs=" << biterrs << " "<< flush;
		 ++blockcount;
		 cout << "blockcount=" << blockcount << " " << flush;
		 cout << "  SNR(dB)=" << SNRdb << " proberr=" << 
			double(biterrs)/double(numbits) << flush;
	  }
	  while(biterrs < numerrstocount);
	  //while(blockcount < nblock);
	  cout << endl << "SNR(dB)=" << SNRdb << " proberr=" << 
		 double(biterrs)/double(numbits) << endl;
   }
}


/*
Local Variables:
compile-command: "g++ -o testturbodec2 -g testturbodec2.cc BCJR.cc BinConvIIR.cc interleave.cc Turboenc.cc Turbodec.cc gran.cc uran.cc -Wno-deprecated"
End:
*/



//  Program: testrepcode.cc
//  Test the repetion code performance
//
//  Todd K. Moon, May 31, 2004
//

#include <iostream>
using namespace std;
#include <math.h>

double gran(void);

int main()
{
   int n = 11;					// number of repetitions in code
   int hardsoft = 1;			// 0 = soft, 1 = hard decoding

   double SNRdBstart = 0;		// set the range of the plot
   double SNRdBend = 10;
   double SNRdBstep = 0.5;
   double SNRdB,SNR;


   double sigma2,sigma;			// noise variance and std. deviation
   int numerrstocount = 100;	// number of errors to count to estimate Pe
   double Ecsqrt = 1;			// coded signal amplitude
   long numbits;				// number of bits total so far
   long biterrs;				// number of bit errors so far
   int i;


   double R = 1./double(n);		// code rate
   int t = (n-1)/2;				// random decoding distance
   int ecount;					// count number of decoding errors

   double *r = new double[n];	// received vector;
   double Recv;					// soft decision statistic

   // loop through the different SNR values
   for(SNRdB = SNRdBstart; SNRdB <= SNRdBend; SNRdB += SNRdBstep) {
	  SNR = pow(10.,SNRdB/10.);	// Convert from dB: SNR=Eb/N0
	  biterrs = 0;				// reset counters for this SNR
	  numbits = 0;
	  sigma2 = Ecsqrt*Ecsqrt/(2.*SNR*R); // compute the variance from the SNR
	  sigma = sqrt(sigma2);

	  // assume that 0 is the input bit (transmit -sqrt(Ec) in each position)

	  // loop until enough detection errors have occured
	  while(biterrs < numerrstocount) {
		 numbits++;				// increment number of bits sent

         // generate the received codevector
		 for(i = 0; i < n; i++) {	
			r[i] = -Ecsqrt + sigma*gran(); // add on noise with variance sigma2
		 }
		 if(hardsoft == 0) {	// do soft decoding
		 // do the soft detection
			Recv = 0;
			for(i = 0; i < n; i++) {
			   Recv += r[i];
			}
			if(Recv > 0) {			// decoded incorrectly
			   biterrs++;
			}
		 }
		 else {					// do hard decoding
			ecount = 0;			// count the number of decoding errors
			for(i = 0; i < n; i++) {
			   if(r[i] > 0) ecount++;
			}
			if(ecount > t) {
			   biterrs++;
			}
		 }
	  }
	  cout << "SNR(dB)=" << SNRdB << " proberr=" << 
		 double(biterrs)/double(numbits) << flush << endl;
   }
}


/*
Local Variables:
compile-command: "g++ -o testrepcode -g testrepcode.cc gran.cc -Wno-deprecated"
End:
*/



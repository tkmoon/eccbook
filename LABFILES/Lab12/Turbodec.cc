// Turbodec.cc -- a Turbo decoder
// Todd K. Moon

#include <math.h>

#include "Turbodec.h"
#include "BPSKmod.h"
#include <iostream>
using namespace std;

Turbodec::Turbodec(int deg, unsigned int h_in, unsigned int g_in, 
				   int in_blocklen, 
				   unsigned int interleaveseed,
				   unsigned char **in_P, int in_puncturelen)
   : enc(deg,h_in,g_in,in_blocklen,interleaveseed),   // build the encoder
     bcjr(enc.Enc1,in_blocklen)
{
   blocklen = in_blocklen;
   numbranch = 2;
   CALLOCMATRIX(prior1,double, blocklen,numbranch);
   CALLOCMATRIX(prior2,double,blocklen,numbranch);
   CALLOCMATRIX(r1,double,blocklen,2);
   CALLOCMATRIX(r2,double,blocklen,2);
   R = 1./3.;					// default rate without puncturing

   // save the puncture information, if there is any
   // Puncture only the parity information
   puncturelen = in_puncturelen;
   if(in_P && puncturelen) {   // if there is puncture information, save the puncture matrix
   // e.g.:  P = [ 1 0
   //              0 1 ]
	  int sumP=0;
	  CALLOCMATRIX(P,unsigned char,2,puncturelen);
	  for(int i = 0; i < 2; i++) {
		 for(int j = 0; j < puncturelen; j++) {
			P[i][j] = in_P[i][j];
			if(P[i][j]) sumP++;
		 }
	  }
	  R = double(enc.Enc1.k*puncturelen)/double(sumP+puncturelen*enc.Enc1.k);
   }
   else {
	  P = 0;
	  puncturelen = 0;
   }
   puncturecycle = 0;  // set the counter that indexes the column puncturing from

}


double**
Turbodec::decode(const double *r,int numit, unsigned int finalstate1,
				 unsigned int finalstate2)
{
   // pull out the different portions of the data
   double *syst = new double [blocklen];
   double *systpi = new double[blocklen];
   int i,j,k;
   int inp;
   double sum,t;
   BPSKmod mod1;
   // Un-commute and/or un-puncture the data
   if(P && puncturelen) {
	  puncturecycle = 0;
	  int datalen = int(blocklen/R + 0.2);
	  for(i = 0,j=0; j < blocklen; j++) {
		 syst[j] = r[i];  // save for de-interleaving
		 r1[j][0] = r[i];
		 i++;
		 if(P[0][puncturecycle]) {  // if Parity 1 used
			r1[j][1] = r[i++];
		 }
		 else {					// else pad with zero
			r1[j][1] = 0;
		 }
		 if(P[1][puncturecycle]) { // if Parity 2 used
			r2[j][1] = r[i++];
		 }
		 else {					// else pad with zero
			r2[j][1] = 0;
		 }
		 puncturecycle = (puncturecycle+1) % puncturelen;
	  }
   }
   else {
	  for(i = 0,j=0; i < 3*blocklen; i+=3,j++) {
		 syst[j] = r[i];
		 r1[j][0] = r[i];
		 r1[j][1] = r[i+1];
		 r2[j][1] = r[i+2];
	  }
   }
   enc.interleaver.Pi(syst,systpi);
   for(i = 0; i < blocklen; i++) {
	  r2[i][0] = systpi[i];
   }

   delete[] syst; delete[] systpi;

   // initialize the probabilities
   double pinit = 1/double(numbranch);
   for(i = 0; i < blocklen; i++) {
	  for(j = 0; j < numbranch; j++) 
		 prior1[i][j] = pinit;
   }

   double **priora, **priorb;	// pointer to probability storage
   double **rptr, **prior;
   cout << "finalstate1=" << finalstate1 << "  finalstate2=" << finalstate2 << endl;
   for(i = 0; i < numit; i++) {

// 	  priora = bcjr.MAP1((const double **)r1,(const double **)prior1,
//  						 finalstate1); // compute extrinsic information

// 	  int b; double t, t1, t2;
// 	  cout << "Map1: " << endl;
// 	  for(j = 0; j < 10; j++) {
// 		 for(b = 0; b < 2; b++) {
// 			cout << priora[j][b] << " ";
// 		 }
// 		 cout << endl;
// 	  }

// 	  priora = bcjr.MAP2((const double **)r1,(const double **)prior1,
//  						 finalstate1); // compute extrinsic information

// 	  cout << "Map2: " << endl;
// 	  for(j = 0; j < 10; j++) {
// 		 for(b = 0; b < 2; b++) {
// 			cout << priora[j][b] << " ";
// 		 }
// 		 cout << endl;
// 	  }
// 	  cout  << "Map2, adjusted" << endl;
// 	  for(j = 0; j < 10; j++) {
// 		 b = 0; t = r1[j][0] - (2.*b-1.);
// 		 t1 = prior1[j][b]*priora[j][b]*exp(-t*t/(2*sigma2));
// 		 cout << "extra fact: " << "sigma2=" << sigma2 << " " <<
// 			prior1[j][b] << " " <<
// 			"r: " << r1[j][0] << " inp: " << b << " " << exp(-t*t/(2*sigma2))
// 			  << " " << prior1[j][b]*exp(-t*t/(2*sigma2)) << endl;
// 		 b = 1; t = r1[j][0] - (2.*b-1.);
// 		 t2 = prior1[j][b]*priora[j][b]*exp(-t*t/(2*sigma2));
// 		 cout << "extra fact: " << prior1[j][b] << " "<< 
// 			"r: " << r1[j][0] << "inpu: " << b << " " << exp(-t*t/(2*sigma2))
// 			  << " " << prior1[j][b]*exp(-t*t/(2*sigma2)) << endl;
// 		 cout << t1/(t1+t2) << " " << t2/(t1+t2) << endl;
// 	  }


	  priora = bcjr.MAP2((const double **)r1,(const double **)prior1,
 						 finalstate1); // compute extrinsic information

	  enc.interleaver.Pi((const double **)priora,prior1,
						 numbranch);
								// permute the extrinsic information
	  if(i==numit-1) break;
	  priorb =  bcjr.MAP2((const double **)r2,(const double **)prior1,
						  finalstate2); 
                                // compute the extrinsic information
	  enc.interleaver.Piinv((const double **)priorb,prior1,
							numbranch);
	  // permute the extrinisic information
   }
   // Include the non-extrinsic information in the final  decoding
   prior = bcjr.MAP1((const double **)r2, (const double **)prior1,finalstate2);
   enc.interleaver.Piinv((const double **) prior,prior1,
						 numbranch);
   prior = prior1;
   return prior;
}

/*
Local Variables:
compile-command: "g++ -c Turbodec.cc"
End:
*/


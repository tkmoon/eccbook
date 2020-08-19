// Turboenc.cc -- a Turbo encoder
// Todd K. Moon

#include "Turboenc.h"
#include <iostream>
using namespace std;

Turboenc::Turboenc(int deg, unsigned int h_in, unsigned int g_in, 
				   int in_blocklen, 
				   unsigned int interleaveseed,
				   unsigned char **in_P, int in_puncturelen)
   : Enc1(1,2,deg,&h_in,g_in),  // build the systematic convolutional coders
     Enc2(1,2,deg,&h_in,g_in),
     interleaver(in_blocklen,interleaveseed)
{
   blocklen = in_blocklen;
   inbits = new unsigned char[blocklen];
   interleavebits = new unsigned char[blocklen];
   R = 1./3.;					// default rate if no puncturing
   // save the puncture information, if there is any
   // The puncturing describes puncturing of parity bits only
   // The systematic part is not punctured
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
	  R = double(Enc1.k*puncturelen)/double(sumP+puncturelen*Enc1.k);
   }
   else {
	  P = 0; 	  puncturelen = 0;
   }
   puncturecycle = 0;  // set the counter that indexes the column puncturing from

   outputbits = new unsigned int[int(blocklen/R)];

}


unsigned int*
Turboenc::encode(const unsigned char *ins)
{
   int i,j,k1;
   unsigned int *codebits;

   Enc1.setstate(0);			// make sure starting from state 0
   Enc2.setstate(0);
   k1 = 0;
   puncturecycle = 0;
// cout <<"States for first encoder:";
   for(i = 0; i < blocklen; i++) {
	  inbits[i] = ins[i];
	  codebits = Enc1.encode(&ins[i]);
//cout << "state1: " << Enc1.getstate() << " ";
	  outputbits[k1++] = codebits[0];
	  if(P) { // if punctured
		 if(P[0][puncturecycle]) { // puncturing on 1st encoder output
			outputbits[k1++] = codebits[1];
		 }
		 if(P[1][puncturecycle]) { // puncturing on 2nd encoder output
			k1++;  // leave space for the parity from the other stream
		 }
		 puncturecycle = (puncturecycle+1) % puncturelen;
	  }
	  else {   // not punctured
//cout << Enc1.getstate() << " ";
		 outputbits[k1++] = codebits[1];  // parity bit
// cout << int(codebits[0]) << int(codebits[1]) << " ";
		 k1++;					    // skip the parity for the other encoder
	  }
   }
 cout << endl;
   interleaver.Pi(inbits,interleavebits);
   if(P) { puncturecycle = 0; k1 = 0;}
   else k1 = 2;
// cout << "inputs/States for second encoder: ";
   for(i = 0; i < blocklen; i++) {
	  codebits = Enc2.encode(&interleavebits[i]);
//cout << "state2: " << Enc2.getstate() << " ";
// cout << "in: " << int(interleavebits[i]) << " state: " << Enc2.getstate() << " ";
	  if(P) { // if punctured
		 k1++;					// skip the systematic bit
		 if(P[0][puncturecycle]) k1++; // skip the parity for the other encoder
		 if(P[1][puncturecycle]) {
			outputbits[k1++] = codebits[1];
		 }
		 puncturecycle = (puncturecycle+1) % puncturelen;
	  }
	  else {
		 outputbits[k1] = codebits[1]; // save the parity in the data stream
//cout << int(codebits[1]) << " ";
		 k1 += 3;
	  }
   }
// cout << endl;
   return outputbits;
}

/*
Local Variables:
compile-command: "g++ -c Turboenc.cc"
End:
*/


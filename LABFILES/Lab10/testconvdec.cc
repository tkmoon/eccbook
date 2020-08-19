//  Program: testconvdec.cc
//
//  Todd K. Moon
//
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include <iostream>
using namespace std;

extern "C" {
#include <stdio.h>
}

#include "BinConvFIR.h"
#include "BinConvdec01.h"
#include "BinConvdecBPSK.h"
#include "matalloc.h"
#include "BPSKmodvec.h"

int main()
{
   int i;
   int decout;

   unsigned int **h;

   int k = 1;					// two inputs
   int n = 2;					// three outputs
   CALLOCMATRIX(h,unsigned int, k,n);
   // G = [D^2+1 D^2+D+1]
   h[0][0] = 5;  h[0][1] = 7;  // first row of G
   int p[] = {2};        // degrees of rows of G

   BinConvFIR conv(k,n,p,h);
   BinConvdec01 decoder(conv,8);

   unsigned char d[] = {1,1,0,0,1,0,1,0};
   int N = sizeof(d)/k;			// number of time steps in test array
   unsigned int *out;
   for(i = 0; i < N; i++) {
	  out = conv.encode(&d[i]);
	  // corrupt the data
	  if(i == 2)
		 out[0] = 0;
	  if(i == 3)
		 out[1] = 0;

	  cout << "in: " << int(d[i]) << " out: ";
	  for(int j = 0; j < n; j++) {
		 cout << int(out[j]) << " ";
	  }
	  cout << " state: " << conv.getstate() << endl;
	  // call the decoder
	  decout = decoder.decode(out);
      // decoder.showpaths();
	  if(decout) {
		 cout << "INPUT FOUND: " << decoder.inputs << endl;
	  }
   }
   cout << endl;
   // dump out the rest
   decoder.showpaths();
   while(decoder.getinpnow(1)) {
	  cout << "INPUT: " << decoder.inputs << endl;
   }


   // test the soft-decision decoder
   BinConvdecBPSK decoderBPSK(conv,8);
   BPSKmodvec modulator(n);
   conv.setstate(0);

   double *mod;					// modulated data
   for(i = 0; i < N; i++) {
	  out = conv.encode(&d[i]);
	  // corrupt the data
	  if(i == 2)
		 out[0] = 0;
	  if(i == 3)
		 out[1] = 0;
	  mod = modulator.mod(out);

	  cout << "in: " << int(d[i]) << " out: ";
	  for(int j = 0; j < n; j++) {
		 cout << int(mod[j]) << " ";
	  }
	  cout << " state: " << conv.getstate() << endl;
	  // call the decoder
	  decout = decoderBPSK.decode(mod);
      // decoderBPSK.showpaths();
	  if(decout) {
		 cout << "INPUT FOUND: " << decoderBPSK.inputs << endl;
	  }
   }
   cout << endl;
   // dump out the rest
   decoderBPSK.showpaths();
   while(decoderBPSK.getinpnow(1)) {
	  cout << "INPUT: " << decoderBPSK.inputs << endl;
   }
}


/*
Local Variables:
compile-command: "g++ -o testconvdec -g testconvdec.cc Convdec.cc BinConvFIR.cc BinConvdec01.cc BinConvdecBPSK.cc"
End:
*/

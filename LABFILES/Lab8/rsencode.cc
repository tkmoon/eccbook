//  Program: rsencode.cc --- encode a data file
//  Todd K. Moon
// 
// rsencode [-t t] infile outfile
// 
// [-t t] -- specify random error correction capability.
//           default = 3

// The data are read with the first byte in the block 
// corresponding to the constant coefficient of m(x).
// In order to handle the last block of the file correctly, 
// the last block is written out with its length.

#include "RSenc.h"
extern "C" {
#include <stdio.h>				// fopen, fread, fwrite functions
#include <string.h>				// strcmp function
#include <stdlib.h>				// atoi function
}

int main(int argc, char *argv[])
{
   GFNUM2m::initgf(8,0x11D);  //   1 0001 1101   x^8 + x^4 + x^3 + x^2 + 1
   GFNUM2m::setouttype(GFvector);
   int i;
   int t = 3;					// default error correction capability
   FILE *fin, *fout;
   if(argc < 3) {
	  cout << "Usage: " << argv[0] << "[-t t] infile outfile\n";
	  exit(-1);
   }
   i = 1;
   if(!strcmp(argv[i],"-t")) { // if t is specified
	  ++i;
	  t = atoi(argv[i++]);
   }
   int n = 255;					// coded block length
   int k = 255-2*t;				// input block length
   RSenc encoder(n,k,t);		// build the encoder object
   unsigned char m[k];			// data space to read in
   unsigned char  c[n];			// encoded data

   if((fin = fopen(argv[i++],"rb")) ==NULL) {// open input file to read
	  cout << "Cannot open file to read\n";
	  exit(-1);
   }
   if((fout = fopen(argv[i],"wb")) == NULL) {  // open output file to write
	  cout << "Cannot open file to write\n";
	  exit(-1);
   }
   int nread, nwritten;
   while((nread = fread(m, 1, k, fin))) {
	  if(nread < k) {  // last block -- not all all k bytes read
		 cout << "last block: nread="<< nread << endl;
		 m[nread] = nread;
		 for(i = nread+1; i < k; i++) {
			m[i] = 0;			// zero out the rest of the block
		 }
	  }
	  encoder.encode(m,c);		// encode the data
	  nwritten = fwrite((void *)c,1,n,fout);
   }
   fclose(fin);
   fclose(fout);
}

/*
Local Variables:
compile-command: "g++ -o rsencode -g rsencode.cc RSenc.cc GFNUM2m.cc"
End:
*/



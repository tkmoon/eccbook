//  Program: rsdecode.cc --- decode a data file
//  Todd K. Moon
// 
// rsencode [-t t] infile outfile
// 
// [-t t] -- specify random error correction capability.
//           default = 3

// In encodding, the data are read with the first byte in the block
// corresponding to the constant coefficient of m(x).  In order to
// handle the last block of the file correctly, the last block is
// written out with its length.  Thus, to decode correctly, the last
// block should use the length and write only the correct number of bytes

#include "RSdec.h"
extern "C" {
#include <stdio.h>				// fopen, fread, fwrite functions
#include <string.h>				// strcmp function
#include <stdlib.h>				// atoi function
}
#include "polynomialT.cc"
// create instantiations of the polynomial class of type GFNUM2m
template class polynomialT<GFNUM2m>;
template class polytemp<GFNUM2m>;

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
   RSdec decoder(t,n);			// build the encoder object
   unsigned char r[n];			// data space to read in
   unsigned char  dec[n];		// decoded data

   if((fin = fopen(argv[i++],"rb")) ==NULL) {// open input file to read
	  cout << "Cannot open file to read\n";
	  exit(-1);
   }
   if((fout = fopen(argv[i],"wb")) == NULL) {  // open output file to write
	  cout << "Cannot open file to write\n";
	  exit(-1);
   }
   int nread, nwritten;
   int prev = 0;
   while((nread = fread(r, 1, n, fin))) {
	  // not the last block --- write out previous message
	  if(prev) {
		 fwrite((void *)(dec+n-k),1,k,fout);
	  }
	  decoder.decode(r,dec);		// decode the data
	  prev = 1;
   }
   // now write out the last block
   // find the last nonzero element
   int blocklen=0;
   for(i = n-1; i >= 0; i--) {
	  if(dec[i]) {
		 if(i==n-1) { // if the block is full
			if(dec[i]==n-1) blocklen = dec[i];
			else blocklen = n-(n-k);
		 }
		 else {
			blocklen = dec[i];
		 }
		 break;
	  }
   }
   cout << "last block: blocklength = " << blocklen << endl;
   fwrite((void *)(dec+n-k),1,blocklen,fout);

   fclose(fin);
   fclose(fout);
}

/*
Local Variables:
compile-command: "g++ -o rsdecode -g rsdecode.cc RSdec.cc GFNUM2m.cc ChienSearch.cc"
End:
*/



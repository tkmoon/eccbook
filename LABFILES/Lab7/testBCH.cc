//  Program: testBCH.cc --- test the BCH decoder
//  Todd K. Moon
#include <iostream>

using namespace std;
#include "BCHdec.h"

// instantiate the berlekamp-massey algorithm for GFNUM coefficients
template <class GFNUM2m> void
berlmass2(const GFNUM2m* s, int n, GFNUM2m* c, int& L);

template <class GFNUM2m> void
berlmass2(const GFNUM2m* s, int n, GFNUM2m* c, int& L);


int main()
{
   int i,j,l;
   int j1,j2,j3;
   GFNUM2m::initgf(4,0x13);  // 1 0011 = 1+d+d^4
   int t = 3;
   int n = 15;
   BCHdec decoder(t,n);
   GF2 r[n];    // the received vector
   GF2 dec[n];    // the decoded vector
   // r[0] = 1; r[1] = 0; r[2] = 0; r[3] = 1; r[4] = 1; r[5] = 0; r[6] = 0;
   // r[7] = 0; r[8] = 0; r[9] = 1; r[10] = 0; r[11] = 0; r[12] = 0; r[13] = 1;
   // r[14] = 0;
   // decoder.decode(r,dec);
   // for(i = 0; i < n; i++) {
   // 	  cout << dec[i] << " ";
   // }
   // exit(-1);
   int allcorrect = 1;
   for(j = 0; j < n; j++) r[j] = 0;  // clear out previous
   for(j1 = 0; j1 < n; j1++) {
	  for(j2 = 0; j2 < n; j2++) {
		 for(j3 = 0; j3 < n; j3++) {
			r[j1] = r[j2] = r[j3] = 1;
			decoder.decode(r,dec);
			for(i = 0; i < n; i++) {
			   if(dec[i] != 0) {
				  allcorrect = 0;
				  cout << "Uncorrected error: (" << j1 <<"," << j2 << ","
					   << j3 << ")" << endl;
				  break;
			   }
			}
			r[j1] = r[j2] = r[j3] = 0;  // clear out previous
		 }
	  }
   }
   if(!allcorrect) {
	  cout << "Errors found!  Decoder not working correctly!" << endl;
   }
   else {
	  cout << "All patterns of 3 errors corrected!" << endl;
   }
}

/*
Local Variables:
compile-command: "g++ -o testBCH -g testBCH.cc BCHdec.cc GFNUM2m.cc ChienSearch.cc"
End:
*/



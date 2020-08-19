//  Program: testRS.cc --- test the RS decoder
//  Todd K. Moon

#include "RSdec.h"
#include "RSenc.h"
#include "polynomialT.cc"

double uran(void);

// create instantiations of the polynomial class of type GFNUM2m
template class polynomialT<GFNUM2m>;
template class polytemp<GFNUM2m>;


int main()
{
   int i,j,l;
   int j1,j2,j3;

   GFNUM2m::initgf(8,0x11D);  //   1 0001 1101   x^8 + x^4 + x^3 + x^2 + 1

   int n = 255;
   int k = 249;
   int t = 3;
   int j0=1;					// starting index
   RSdec decoderb(t,n);
   GFNUM2m r[255];
   GFNUM2m decb[255];
   int loc[3];
   int val[3];
   int nerror;
   int allcorrect = 1;

   // test 100 received vectors
   for(i = 0; i < 100; i++) {
	  cout << "i=" << i << endl;

	  // generate a vector with a random number of errors
	  nerror = int(3*uran()) + 1;  // number of errors 1 or 2 or 3
	  for(j = 0; j < nerror; j++) {   
		 loc[j] = int(n*uran());  // generate error locations
		 val[j] = int(256*uran()); // generator error values
		 r[loc[j]] = val[j];
	  }

	  // Call the decoder here...
	  // (fill in the blanks)

	  // make sure all decoded values are correct
	  for(j = 0; j < n; j++) {
		 if(decb[j] != 0) {
			allcorrect = 0;
			cout << "Undecoded error: i=" << i << " " << loc[0] << " " <<
			   loc[1] << " " << loc[2] << " " << val[0] << " " << val[1] <<
			   " " << val[2] << endl;
			break;
		 }
	  }
	  for(j = 0; j < nerror; j++) {
		 r[loc[j]] = 0;
	  }
   }
   if(!allcorrect) {
	  cout << "not all errors correct!" << endl;
   }
   else {
	  cout << "All errors correted!" << endl;
   }
}

/*
Local Variables:
compile-command: "g++ -o testRS -g testRS.cc GFNUM2m.cc RSdec.cc ChienSearch.cc RSenc.cc uran.cc"
End:
*/



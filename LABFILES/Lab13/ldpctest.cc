// ldpctest.cc -- test the low-density parity-check code
// decoder


#include "ldpcdec.h"
#include <iostream>
#include <math.h>
using namespace std;

int main()
{
   int i;
   int numloops;
   int success;
   double a = 2;
   double sigma2 = 2;
   LDPCDEC ldpcdec("Asmall.txt",1,1);
   ldpcdec.setsigma2(sigma2);
   ldpcdec.setsigamp(a);
   ldpcdec.printsparseA();
   double p1[]  = {.22, .16, .19, .48,.55, .87,.18,.79,.25,.76}; 
   double y[10];
   for(i = 0; i < 10; i++) {
	  y[i] = log(1./p1[i]-1)/(-2*a)*sigma2;
   }
   VECDUMP(p1,ldpcdec.N);
   VECDUMP(y,ldpcdec.N);
   numloops = 0;
   success = ldpcdec.decode(y,1,10,numloops);  // decode and print
   cout << "Success=" << success << 
	  "  number of loops=" << numloops << endl;
}

/*
Local Variables:
compile-command:"g++ -o ldpctest -g ldpctest.cc  ldpcdec.cc"
End:
*/

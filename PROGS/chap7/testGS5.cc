//  Program: testGS.cc --- test the GS decoder and multivar poly arithmetic
//  This program actually goes through a decoding process for a RS code

//  Todd K. Moon, Feb 7, 2004


#include <math.h>
#define TYPE ModAr
// This gives a sort of cheap template capability, so that 
// examples can be tested in either GF(2^m) or GF(p)

#include "ModAr.h"

#include "polynomialT.cc"


// create instantiations of the polynomial class 
template class polynomialT<TYPE>;
template class polytemp<TYPE>;
template class polynomialT<polynomialT<TYPE> >;

polynomialT<TYPE>
kotter1(int n,int k,TYPE *ai, TYPE *bi,int *wdeg,int rlex);

int main()
{
   TYPE::setdefaultm(5);		// work over GF19
   TYPE::setshowmod(0);		// don't show modulus

   int n = 5;
   int k = 2;
   int d = n-k+1;
   int t = int(floor((d-1)/2));
   
   cout << "n=" << n << "  k=" << k << "  d=" << d << "  t0=" << t << endl;

   POLYC(TYPE,m,{1,4});
   cout << "message polynomial=" << m << endl;
   TYPE cc[n];				// code vector
   TYPE sc[n];				// support set
   int supset1 = 0, supset2 = 4; // support set of code : 0..4
   int i,j;
   for(j=0, i = supset1; i <= supset2; i++, j++) {
	  cc[j] = m(i);
	  sc[j] = i;
   }
   polynomialT<TYPE> c(n-1,cc);
   polynomialT<TYPE> ss(n-1,sc);
   cout << "code polynomial=" << c << endl;
   // build error 
   int nerr = 1;
   int errloc[n];
   TYPE errval[n];
   errloc[0] = 1;
   errloc[1] = 3;
   errloc[2] = 5;
   errloc[3] = 7;
   errloc[4] = 9;

   errval[0] = 2;
   errval[1] = 2;
   errval[2] = 2;
   errval[3] = 2;
   errval[4] = 2;
   TYPE ec[n];
   for(i = 0; i < nerr; i++) {
	  ec[ errloc[i]] = errval[i];
   }
   polynomialT<TYPE> e(n-1,ec);
   cout << "e=" << e << endl;
   polynomialT<TYPE> r(c+e);	
   cout << "received polynomial=" << r << endl;
   int wdeg[2] = {1,k-1};		// weighted degree
   int rlex = 1;				// rlex=1 for rlex; rlex=0 for lex


   polynomialT<TYPE> p;			// this will be the message poly
   p = kotter1(n,k,ss.getarray(),r.getarray(),wdeg,rlex);
   cout << "p=" << p << endl;
}



#include "kotter1.cc"  

/*
Local Variables:
compile-command: "g++ -o testGS5 -Wno-deprecated -g testGS5.cc ModAr.cc polynomialT.cc"
End:
*/



//  Program: testGS.cc --- test the GS decoder and multivar poly arithmetic

//  Todd K. Moon, Jan 24, 2004

#define TYPE GFNUM2m

#include "GFNUM2m.h"
#include "polynomialT.cc"

#define BIGINT 65535

// create instantiations of the polynomial class of type GFNUM2m
template class polynomialT<GFNUM2m>;
template class polytemp<GFNUM2m>;

template class polynomialT<polynomialT<GFNUM2m> >;

GFNUM2m evaluate(polynomialT<polynomialT<GFNUM2m> > &Q,GFNUM2m a,
				 GFNUM2m b);
polynomialT<polynomialT<GFNUM2m> >
kotter(int n,int L,GFNUM2m *ai, GFNUM2m *bi,int *mi, int *wdeg,int rlex,
	   int m1);


int main()
{
   GFNUM2m::initgf(4,0x13);  // 1 0011 = 1+d+d^4

   // Set up data for the Kotter algorithm test
   int n = 5;					// code length
   int mi[] = {1,1,1,1,1};		// orders at each pt.
   GFNUM2m ai[] = {A^0,A,A^2,A^3,A^4};
   GFNUM2m bi[] = {A^3,A^4,A^5,A^7,A^8};
   int L = 4;
   int wdeg[2] = {1,2};			// weighting for monomial order
   int rlex;					// rlex=1 for rlex; rlex=0 for lex

   // Kotter algorithm
   polynomialT<polynomialT<GFNUM2m> > P;
   P = kotter(n,L,ai,bi,mi, wdeg,rlex,0);
   cout << "interpolating polynomial=" << P << endl;
}

#include "kotter.cc"

/*
Local Variables:
compile-command: "g++ -o testGS1 -Wno-deprecated -g testGS1.cc GFNUM2m.cc polynomialT.cc"
End:
*/



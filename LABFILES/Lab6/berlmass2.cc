//  Program: berlmass2.cc --- implement the Berlekamp-Massey algorithm
//  using arrays, not polynomials

#include "polynomialT.h"

#define MAX2(a,b) (a> b? a:b)

template <class T> void
berlmass2(const T* s, int n, T* c, int& L)
// s = input coefficients s[0],s[1],... s[n-1]
// c = connection polynomial coefficients.  Must be allocated prior to calling
// L = degree of connection polynomial
{
   L = 0;
   T t[n];						// temporary
   T p[n];						// previous
   c[0] = T(1);
   p[0] = T(1);
   int shift = 1;				// n-m
   T dm = T(1);					// previous discrepancy
   T d;							// current discrepancy
   T ddm;						// d/dm;
   int j;
   int oldL;
   for(int k=0; k < n; k++) {

	  // fill in the blanks
	  // ....

	  
      cout << "k=" << k << "  S=" << s[k] << "  d=" << d << "  c=";
	  for(j = 0; j<= L; j++) cout << c[j] << " ";
	  cout << "  L=" << L << "  p="; for(j=0;j<=oldL;j++) cout << p[j] <<" ";
	  cout << " l=" << shift << "  dm=" << dm << endl;
	  cout << "k=" << k << "  d=" << d << "  c=";
	  for(j = 0; j<= L; j++) cout << c[j] << " "; cout << "  L=" << L << endl;

   }
}

#undef MAX2
/*
Local Variables:
compile-command: "g++ -g berlmass.cc"
End:
*/



//  Program: berlmass.cc --- implement the Berlekamp-Massey algorithm
//  Todd K. Moon

#ifndef BERLMASS_CC
#define BERLMASS_CC

#include "polynomialT.h"

template <class T> polynomialT<T>
berlmass(const T* s, int n)
// s = input coefficients s[0],s[1],... s[n-1]
// returns = final connection polynomial
{
   int L = 0;
   polynomialT<T> c;
   c[0] = T(1);
   // cout << "c=" << c << endl;
   polynomialT<T> p(c), t;      // p is previous, t is temporary
   int shift = 1;               // n-m
   T dm = T(1);                 // previous discrepancy
   T d;                         // current discrepancy
   int j;
   for(int k=0; k < n; k++) {

	  // fill in the blanks ....
	  
cout<<"k="<<k << " s=" << s[k] << " d=" << d << "  c=" << c << "  L=" << L;
cout << " p=" << p << "  l=" << shift << " dm=" << dm << endl;

   }
   return c;
}


#endif
/*
Local Variables:
compile-command: "g++ -g berlmass.cc"
End:
*/



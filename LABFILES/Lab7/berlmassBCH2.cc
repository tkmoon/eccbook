//  Program: berlmassBCH2.cc --- implement the Berlekamp-Massey algorithm
//  using arrays, not polynomials,
//  and using the speedups for BCH-type syndromes
//  Todd K. Moon

#include "polynomialT.h"

#define MAX2(a,b) (a> b? a:b)

template <class T> void
berlmassBCH2(const T* s, int n, T* c, int& L)
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
   for(int k=0; k < n; k+=2) {
	  d = s[k];
	  for(j=1; j <= L; j++) {
		 d += c[j]*s[k-j];
	  }
	  if(d == 0) {
		 shift++;				// increment time since last shift
	  }
	  else {
		 ddm = d/dm;
		 if(2*L > k) {			// update without length change
			//c = c - (p << shift)*(d/dm);
			for(j = shift; j <= L; j++) {
			   c[j] -= p[j-shift]*ddm;
			}
			shift++;
		 }
		 else {					// update with length change
			for(j = 0; j <= L; j++) t[j] = c[j];  // t = c;
			oldL = L;
			L = k+1-L;			// new degree
			for(j = shift; j <= oldL; j++) {  // c = c - (p<< shift)*(d/dm);
			   c[j] -= p[j-shift]*ddm;
			}
			for(j = oldL+1; j < shift; j++) c[j] = 0;
			for(j = MAX2(shift,oldL+1); j <= L; j++) {  // c = c - (p<< shift)*(d/dm);
			   c[j] = -p[j-shift]*ddm;
			}
			for(j = 0; j <= oldL; j++) { p[j] = t[j]; } // p = t;
			dm = d;
			shift = 1;
		 }
	  }
	  shift++; // account for odd k skipped
	  //      cout << "k=" << k << "  S=" << s[k] << "  d=" << d << "  c=";
	  // for(j = 0; j<= L; j++) cout << c[j] << " ";
	  // cout << "  L=" << L << "  p="; for(j=0;j<=oldL;j++) cout << p[j] <<" ";
	  // cout << " l=" << shift << "  dm=" << dm << endl;
   }
}

#undef MAX2
/*
Local Variables:
compile-command: "g++ -g -c berlmassBCH2.cc"
End:
*/



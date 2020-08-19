// testpoly1.cc -- test the polynomialT<TYPE> class using
// double and ModAr types

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "polynomialT.cc"
#include "ModArnew.h"

#include <iostream>
using namespace std;

// make some typedefs for shorthand
typedef ModAr<5> Z5;
typedef ModAr<2> Z2;
typedef polynomialT<double> polyR;
typedef polynomialT<Z5> polyZ5;
typedef polynomialT<Z2> polyZ2;

// cause template classes to instantiate with the desired types
template class polynomialT<double>;
template class polynomialT<Z5>;
template class polynomialT<Z2>;

int main()
{
   // First check polynomials in R[x]

   // initialize some polynomials

   cout << "Polynomials with double coefficients\n";

#ifdef __GNUC__   
   POLYC(double,p1,{1,0,0,1,3,0,4,3}); // 1+x^3 + 3x^4 + 4x^6 + 3x^7
   POLYC(double,p2,{0,1,0,1,3});       // x+x^3+ 3x^4
#else  // for microsoft compilers
   polyR p1; { double p1a[] = {1,0,0,1,3,0,4,3}; p1.setc(p1a,7);  }
   polyR p2; { double p2a[] = {0,1,0,1,3}; p2.setc(p2a,4);  }
#endif

   cout << "p1=" << p1 << endl;
   cout << "p2=" << p2 << endl;

   polyR q = p1/p2;
   polyR r = polyR::getlastremainder();
   cout << "quotient (p1/p2): " << q << endl;
   cout << "remainder (p1/p2): " << r << endl;
   cout << "sum: p1+p2=" << p1+p2 << endl;
   cout << "difference: p1-p2=" << p1-p2 << endl;
   cout << "product: p1*p2=" << p1*p2 << endl;
   polyR s = p1/p2 + p1*(p2+p2)*p1*p1;
   cout << "s = p1/p2 + p1*(p2+p2)*p1*p1=" << s << "\n";

   // Now check polynomials in Z_5[x]


   cout << "Polynomials with Z5 coefficients\n";
#ifdef __GNUC__
   POLYC(Z5,p1m,{1,0,0,1,3,0,4,3}); // 1+x^3 + 3x^4 + 4x^6 + 3x^7
   POLYC(Z5,p2m,{0,1,0,1,4});       // x+x^3+4x^4
#else  // for microsoft compilers
   polyZ5 p1m; { Z5 p1a[] = {1,0,0,1,3,0,4,3}; p1m.setc(p1a,7);  }
   polyZ5 p2m; { Z5 p2a[] = {0,1,0,1,4}; p2m.setc(p2a,4);  }
#endif
   cout << "p1m=" << p1m << endl;
   cout << "p2m=" << p2m << endl;
   polyZ5 qm = p1m/p2m;
   polyZ5 rm = polyZ5::getlastremainder();
   cout << "quotient (p1/p2): " << qm << endl;
   cout << "remainder (p1/p2): " << rm << endl;
   cout << "sum: p1+p2=" << p1m+p2m << endl;
   cout << "difference: p1-p2=" << p1m-p2m << endl;
   cout << "product: p1*p2=" << p1m*p2m << endl;
   polyZ5 sm = p1m/p2m + p1m*(p2m+p2m)*p1m*p1m;
   cout << "s = p1/p2 + p1*(p2+p2)*p1*p1=" << sm << "\n";

   // Now a polynomial in GF(2);

   cout << "Polynomials with Z2 coefficients\n";
#ifdef __GNUC__
   POLYC(Z2,x,{0,1});		// x
   POLYC(Z2,den,{1,1,0,1,1,0,1}); // 1 + x + x^3 + x^4 + x^6
#else
   polyZ5 x; { Z2 xa[] = {0,1}; x.setc(xa,1);  }
   polyZ5 den; { Z2 dena[] = {1,1,0,1,1,0,1}; den.setc(dena,6);  }
#endif

   polyZ2 num = (x<<62) + Z2(1);
   cout << "num=" << num << endl;
   cout << "den=" << den << endl;
   cout << "num/den=" << (num/den) << endl;
   cout << "num%den=" << num % den << endl;
   cout << "den/num=" << den/num << endl;
   cout << "den % num=" << den % num << endl;
}

/*
Local Variables:
compile-command:"g++ -o testpoly1 -g testpoly1.cc"
End:
*/

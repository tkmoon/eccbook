//
//  Program: testpolygcd.cc

#include <math.h>
#include "ModArnew.h"
#include "polynomialT.cc"
#include "gcdpoly.cc"

// create instantiations of the polynomial class of type ModAr
template class polynomialT<ModAr<5> >;
template class polytemp<ModAr<5> >;
// create instantiations of the polynomial class of type double
template class polynomialT<double>;
template class polytemp<double>;

typedef ModAr<5> Z5;

// Create an instantiation of the gcd function, for polynomials of type ModAr
template <class Z5 > void
gcd(const polynomialT<Z5> &a, const polynomialT<Z5> &b, 
	polynomialT<Z5> &g,
	polynomialT<Z5> &s, polynomialT<Z5> &t, int sdeg);

// Create an instantiation of the gcd function, for polynomials of type double
// template <double> void
// gcd(const polynomialT<double> &a, const polynomialT<double> &b, 
// 	polynomialT<double> &g,
// 	polynomialT<double> &s, polynomialT<double> &t, int sdeg);


int main()
{
   POLYC(Z5,a,{1,0,0,1,3,0,4,3}); // 1+x^3 + 3x^4 + 4x^6 + 3x^7
   POLYC(Z5,b,{0,1,0,1,4});       // x+x^3+x^4

   polynomialT<Z5> s, t, g;
   a.setprintdir(1);

   cout << "a=" << a << endl;
   cout << "b=" << b << endl;
   gcd(a,b,g,s,t);
   cout << "g=" << g << endl; 
   cout << "s=" << s << endl;
   cout << "t=" << t << endl;

   cout << "\n\n\n";

   double d1d[4] = {2,8,10,4};
   double d2d[4] = {1,7,14,8};

   POLYC(double,ad,{2,8,10,4});  // 2 + 8x + 10x^2 + 4x^3
   POLYC(double,bd,{1,7,14,8});  // 1 + 7x + 14x^2 + 8x^3
   polynomialT<double> sd,td,gd;
   cout << "a=" << ad << endl;
   cout << "b=" << bd << endl;
   gcd(ad,bd,gd,sd,td,0);
   cout << "g=" << gd << endl; 
   cout << "s=" << sd << endl;
   cout << "t=" << td << endl;
}

/*
Local Variables:
compile-command:"g++ -o testpolygcd -g testpolygcd.cc"
End:
*/



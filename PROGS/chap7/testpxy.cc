//  Program: testpxy.cc -- Demonstrate some concepts relating to 
// polynomials with two variables.

//  Todd K. Moon, Feb 12, 2004

// This short program illustrates some concepts relating to 
// using the polynomialT class to implement polynomials with two variables.
// A polynomial Q(x,y) in F[x,y] can be thought of as
// being in (F[x])[y], that is, as polynomials in the variable y
// with coefficients in F[x].  
// The polynomialT class, with its templatized coefficients, can be 
// used to implement this concept.
// For example, the declaration
//
// polynomialT<polynomialT<GFNUM2m> > Q
//
// creates a polynomial whose coefficients are, themselves,
// polynomials, with coefficients in GFNUM2m.
//
// For purposes of interpreting and printing these polynomials, it 
// is helpful to distinguish between the "outer" variable (the y variable
// in the example above, and the inner variable (the x variable).
// The member function setvarname sets the variable name that is printed. 
// For example, 
// 
// Q.setvarname("y")
// 
// sets the variable name of the "outer" variable to y.
// It is also helpful for legibility to have the coefficients (that is, the 
// "inner" polynomials in the variable x, printed with parentheses
// around them.  
// The member function setbeforeafterprint("y","(",")");
// indicates that when printing the _coefficients_ of terms in the 
// variable y to preceed them with "(" and follow them with ")":
//
// Q.setbeforeafterprint("y","(",")");
//
//
// You can think of either the inner or the outer variables as being
// x or y.  For the purposes to which this has been put (The 
// Kotter algorithm and the Roth-Ruckenstein algorithm in conjunction
// with the Guruswami-Sudan algorithm, it is convenient to think of
// the outer variable as y.
//
// Because of the way that the templates work, polynomial operations
// on multi-variable polynomials work just as expected.  That is, these
// can be added, multiplied, etc.
// However, the name conventions established using setvarname and
// setbeforeafterprint are not propagated through the 
// arithmetic operations (for purposes of speed).  It may be
// necessary to set them before display.

#include "GFNUM2m.h"
#include "polynomialT.cc"

// create instantiations of the polynomial class 
template class polynomialT<GFNUM2m>;
template class polytemp<GFNUM2m>;
template class polynomialT<polynomialT<GFNUM2m> >;


main()
{
   GFNUM2m::initgf(4,0x13);  // 1 0011 = 1+d+d^4

   // build some polynomials as coefficients "coefficients" 
   POLYC(GFNUM2m,mc2,{A,A^2});
   cout << "mc2=" << mc2 << endl;
   POLYC(GFNUM2m,mc3,{A^3,A^4});
   cout << "mc3=" << mc3 << endl;
   POLYC(GFNUM2m,mc4,{A^5,A^6});
   cout << "mc4=" << mc4 << endl;

   // Build a polynomial that has these polynomials as
   // coefficients.
   POLYC(polynomialT<GFNUM2m >,Q,{mc2,mc3,mc4});
   Q.setvarname("y");			// set the way Q is printed
   Q.setbeforeafterprint("y","(",")");
   cout << "Q=" << Q << endl;

   // Build some more polynomials
   POLYC(GFNUM2m,p1,{A^5,A^6,A^7});
   POLYC(GFNUM2m,p2,{A^7,A^8,A^9});
   POLYC(GFNUM2m,p3,{A^10,A^11,A^12});

   // Build a polynomial with these polynomial coefficients
   polynomialT<polynomialT<GFNUM2m > > P; // build this one the direct way
   P.setc(0,p1);
   P.setc(1,p2);
   P.setc(2,p3);
   P.setvarname("y");			// set the way P is printed
   P.setbeforeafterprint("y","(",")");
   cout << "P=" << P << endl;

   // Now do some arithmetic
   cout << "Q+P=" << Q+P << endl;
   // This does not preserve the printing conventions, so it is 
   // a little hard to read.
   // So we make a polynomial to hold the result
   polynomialT<polynomialT<GFNUM2m> > result;
   result = Q+P;
   result.setvarname("y");			// and set the way result is printed
   result.setbeforeafterprint("y","(",")");
   cout << "Q+P=" << result << endl;
   result = Q*P;
   result.setvarname("y");			// must reset the printing
   result.setbeforeafterprint("y","(",")");
   cout << "Q*P=" << result << endl;

}

/*
Local Variables:
compile-command: "g++ -o testpxy -Wno-deprecated -g testpxy.cc polynomialT.cc GFNUM2m.o "
End:
*/



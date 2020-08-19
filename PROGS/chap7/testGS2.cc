//  Program: testGS2.cc --- test the GS decoder and multivar poly arithmetic

//  Todd K. Moon, Jan 24, 2004

// the TYPE is defined as a sort of cheap templatized trick,
// so that examples in both GF(2^m) and GF(p) can be setup,
// without the programming overhead of a templatized function.


// Select one of the following two lines
// #define GFTYPE
#define MODARTYPE


#ifdef MODARTYPE
#include "ModAr.h"
#define TYPE ModAr
#endif
#ifdef GFTYPE
#include "GFNUM2m.h"
#define TYPE GFNUM2m
#endif

using namespace std;

#include "polynomialT.cc"
#include "rothruck.h"

// create instantiations of the polynomial class 
template class polynomialT<TYPE>;
template class polytemp<TYPE>;
template class polynomialT<polynomialT<TYPE> >;  // two-variable polynomials

int main()
{
#ifdef MODARTYPE
   TYPE::setdefaultm(5);		// work over GF19
   TYPE::setshowmod(0);		// don't show modulus
#endif
#ifdef GFTYPE
   GFNUM2m::initgf(4,0x13);  // 1 0011 = 1+d+d^4
#endif

   // build polynomial Q(x,y) with given y-roots
   TYPE p1c[] = {1,2};
   polynomialT<TYPE> p1(1,p1c);
   TYPE p2c[] = {4,3,2};
   polynomialT<TYPE> p2(2,p2c);
   TYPE p3c[] = {1,3,1,2};
   polynomialT<TYPE> p3(3,p3c);
   TYPE fc[] = {4,2};
   polynomialT<TYPE> f(1,fc);
   polynomialT<polynomialT<TYPE> > Qnew(polynomialT<TYPE>(1));
   polynomialT<polynomialT<TYPE> > l1;
   l1.setc(1,polynomialT<TYPE>(1));
   l1[0] = -p1;
   Qnew *= l1;
   l1[0] = -p2;
   Qnew *= l1;
   l1[0] = -p3;
   Qnew *= l1;
   Qnew *= f;
   Qnew.setvarname("y");
   Qnew.setbeforeafterprint("y","(",")");
   cout << "Qnew=" << Qnew << endl;
   rpolynode *rptr,*rptr1;		// pointer to list of polynomials
   int D = 2;
   rptr = rothruck(Qnew,D);			// find the y-roots of Q to to degree D
   cout << "y-roots of Q(x,y): " << endl;
   for(rptr1 = rptr; rptr1 != NULL; rptr1 = rptr1->next) {
	  cout << rptr1->f << endl;
   }
}

#include "rothruck.cc"


/*
Local Variables:
compile-command: "g++ -o testGS2 -Wno-deprecated -g testGS2.cc ModAr.cc GFNUM2m.cc polynomialT.cc"
End:
*/



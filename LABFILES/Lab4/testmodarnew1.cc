//******************************************************************
// testmodar1.cc -- test the modularized arithmetic
// Created by Todd K. Moon, Electrical and Computer Engineering Dept.
// Utah State University.    Copyright January 1999
// *****************************************************************

#include "ModArnew.h"

#include <iostream>
using namespace std;

int main()
{
   // ModAr::setdefaultm(5);		// set the default modulo
   // ModAr::setshowmod(0);		// set so modulo is not displayed

   ModAr<5> a(3);
   ModAr<5> b(4);
   ModAr<5> c;
   int aint;
   c = a+b;
   cout << "a=" << a << "  b=" << b << "  c=a+b=" << c << "\n";
   a = 2;
   cout << "a=" << a << "\n";
   c = a+1;
   cout << "a+1=" << a+1 << "\n";
   cout << "1+a=" << 1+a << "\n";
   cout << "c=" << c << "\n";
   aint = a.toint();
   ModAr<5> f = -a;
   cout << "f=-a= " << f << "\n";
   cout << "f^5= " << (f^5) << "\n";
   cout << "a-b= " << a-b << "\n";
   cout << "3-b= " << 3-b << "\n";
   cout << "a*b= " << a*b << "\n";
   cout << "a/b= " << a/b << "\n";
   cout << "1/a= " << 1/a << "\n";
   cout << "3/a= " << 3/a << "\n";
   cout << "a==2 " << (a==2) << "\n";
   cout << "a!=2 " << (a!=2) << "\n";
   cout << "2==a " << (2==a) << "\n";
   cout << "2!=a " << (2!=a) << "\n";
   cout << "a+=b " << (a += b) << "\n";
   cout << "a-=b " << (a -= b) << "\n";
   cout << "a=" << a << " b=" << b << "\n";
   cout << "a*=b " << (a *= b) << "\n";
   cout << "a/=b " << (a /= b) << "\n";
   
   // Now do a few tests mod 6
   ModAr<6> h(3), i(4), j(5);
   cout << "h*i=" << h*i << "\n";
   cout << "h/i=" << h/i << "\n";  // this should be an error and abort
   cout << "all done\n";

}

/*
Local Variables:
compile-command:"g++ -o testmodarnew1 testmodarnew1.cc"
End:
*/

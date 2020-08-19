//  Program: testgfnum.cc
//  test basic operations on GFNUM2m class

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "GFNUM2m.h"

int main()
{
   int g = 0x13;  // 1 0011 = 1+x+x^4
   GFNUM2m::initgf(4,g);
   GFNUM2m a,b,c;
   GF2 a2,b2,c2;
   // GFNUM2m::setouttype(vector);
   a = A^3;
   cout << a << endl;
   b = A^5;
   cout << b << endl;
   cout << "a+b=" << a+b << endl;
   cout << "a+=b=" << (a+= b) << endl;
   a = A^3;
   cout << "a+A^9=" << a+10 << endl;
   cout << "a+A^9=" << 10+a << endl;
   cout << "a=" << a << "  b=" << b << endl;
   cout << "a*b=" << a*b << endl;
   cout << "a*=b=" << (a*=b) << endl;
   a = A^3;
   cout << "a*8=" << (a*8) << endl;
   cout << "a*8=" << (8*a) << endl;
   cout << "a*=8=" << (a*=8) << endl;
   a = A^3;
   cout << "a*A^5=" << a*(A^5) << endl;
   cout << "a^4=" << (a^4) << endl;
   cout << "a^{-4}=" << (a^(-4)) << endl;
   cout << "a^=4=" << (a^=4) << endl;
   a = A^3;
   cout << "a^=(-4)=" << (a^=(-4)) << endl;
   a = A^3;
   cout << "a/b=" << a/b << endl;
   cout << "a/=b=" << (a/=b) << endl;
   a = A^3;
   cout << "a/6=" << a/6 << endl;
   cout << "a/=6=" << (a/=6) << endl;
   a = A^3;
   cout << "1/a=" <<  (1/a) << endl;
   a = A^3;
   cout << "a==b" << (a==b) << endl;
   cout << "a!=b" << (a!=b) << endl;
   cout << "-a=" << -a << endl;
   a2 = 1;
   b = a+a2;
   cout << "a2=" << a2 << "  b=" << b << endl;
}

/*
Local Variables:
compile-command: "g++ -g -o testgfnum testgfnum.cc GFNUM2m.cc"
End:
*/



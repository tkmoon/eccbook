// ******************************************************************
// ModAr.cc -- Modulo Arithmetic
// Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "ModAr.h"
#include <stdlib.h>
#include <math.h>

// Define and set the static variables
int ModAr::defaultm = 2;  // default modulus
int ModAr::showmod = 1;   // show modulus when printing

ModAr 
ModAr::operator^(const int e) const
{
   // Returns a^e.  Caution: the precedence of this operator does not follow
   // conventional mathematics.  Must use with parenthesis, as in a*(b^3).
   int ehere = e;
   int p=1;
   // checking if e==0 is faster than checking i<=e
   // if p==0, there is no need to go further
   for(; ehere&&p; ehere--)
	  p = (p*v) % m;  // mod at each step just to keep growth in check
   // (this is slower than it needs to be)
   return ModAr(p,m);
}

ModAr 
ModAr::operator/(const ModAr &b)  const
{
   // assert(m == b.m);
   // mod operation is done in constructor
   if(b.v==0) {
	  cout << "Division by zero not allowed\n";
	  exit(-1);
   }
   int x,y,g;
   gcd(m,b.v,x,y,g);
   if(g != 1) {
	  cerr << "Error: Attempting to divide by noninvertible number\n";
	  exit(-1);
   }

   return ModAr(v*y,m);
}

ModAr&
ModAr:: operator/=(const ModAr &a) 
{
   // assert(m == a.m);
   if(a.v==0) {
	  cout << "Division by zero not allowed\n";
	  exit(-1);
   }
   int x,y,g;
   gcd(m,a.v,x,y,g);
   if(g != 1) {
	  cerr << "Error: Attempting to divide by noninvertible number\n";
	  exit(-1);
   }
   v = (v * y) % m;
   if(v < 0) v += m;
   return *this;
}

void ModAr::gcd(const int m,const int b,int &x,int &y,int &g) const
{
   // Compute the inverse of b by finding the gcd(m,b)
   int u[3],v1[3],t[3];
   int *up, *vp, *tp, *temp;
   int q;
   u[0] = 1;  u[1] = 0;  u[2] = m;
   v1[0] = 0;  v1[1] = 1;  v1[2] = b;

   up = u;  vp = v1;  tp = t;
   while(vp[2]) {
	  q = up[2]/vp[2];
	  tp[0] = up[0]-vp[0]*q;   // t = u - v*q
	  tp[1] = up[1]-vp[1]*q;
	  tp[2] = up[2]-vp[2]*q;
	  temp = up;
	  up = vp;                // u = v;
	  vp = tp;                // v = t;
	  tp = temp;
   }
   x = up[0];  y=up[1];  g=up[2];
}

ostream&
operator<<(ostream &os, const ModAr &a) 
{
   if(ModAr::showmod) 
	  return os << a.v << '(' << a.m << ')';
   return os << a.v;
}


/*
Local Variables:
compile-command: "g++ -c -g ModAr.cc"
End:
*/

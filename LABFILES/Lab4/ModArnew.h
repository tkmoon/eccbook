// ******************************************************************
// MODAR.h -- MODulo ARithmetic
// Created by Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

// This version of the code uses a template for the modulus.
// This eliminates the need for some run-time checking and storage.
// the compiler can verify some things.
 

#ifndef ModArnew_H
#define ModArnew_H

#include <iostream>
using namespace std;

template <int M>
class ModAr {
protected:
   int v;            // the modulo and value
   static int showmod;	// whether to display the modulo value
   static void gcd(const int m,const int b,int &x,int &y,int &g);
           // Static member function used to compute inverses
public:
	// default constructor
   ModAr(int inv = 0) {
	  v = inv % M;  // Note: m should be >= 2: perhaps add test
	  if(v < 0) 
		 v += M;
   }
	// copy constructor
   ModAr(const ModAr& in) { v = in.v; }

	// destructor
   ~ModAr() {};	

   int getv() const { return v; };
   int getm() const { return M; };
   int character() const { return M; };
	// the int cast operator isn't included because it is identical
	// to the getv function
   int toint() { return (int)v; };  

   static void setshowmod(int inm) { showmod = inm; };
   static int getshowmod() { return showmod; };
   
   // ***********************************************************************
   //                   the overloaded operators
   // ***********************************************************************

   ModAr& operator+=(const ModAr &a) {
	  v = (v + a.v) % M;
	  if(v < 0)  v += M;
	  return *this;
   }
   ModAr& operator-=(const ModAr &a) {
	  v = (v - a.v) % M;
	  if(v < 0)  v += M;
	  return *this;
   }
   ModAr& operator*=(const ModAr &a) {
	  v = (v * a.v) % M;
	  return *this;
   }
   ModAr& operator/=(const ModAr &a);  // defined in ModAr.cc

	// unary -
   ModAr operator-() const { return ModAr(M - v); }

   ModAr operator+(const ModAr &b) const {return ModAr(v+b.v);}
   ModAr operator-(const ModAr &b) const {return ModAr(v-b.v);}
   ModAr operator*(const ModAr &b) const {return ModAr(v*b.v);}
   ModAr operator/(const ModAr &b) const; // defined in ModAr.cc

	// exponentiation operator: watch the precedence!!
	// const int e is changed to int e because there is no
	// way the variable location the user supplies for e is
	// going to be modified.  The function is more efficient
	// if it can modify its own copy of e.
   ModAr operator^(const int e) const; // defined in ModAr.cc

   int operator==(const ModAr &b) const {return (v == b.v);}
   int operator!=(const ModAr &b) const {return (v != b.v);}

   //friend ostream& operator<<(ostream &os, const ModAr<M> &a);
};

template<int M> inline ModAr<M> 
operator+(const int a, const ModAr<M> b) { return b+a;}

template<int M> inline ModAr<M> 
operator-(const int a, const ModAr<M> b) { return ModAr<M>(a)-b;}

template<int M> inline ModAr<M> 
operator*(const int a, const ModAr<M> b) { return b*a;}

template<int M> inline ModAr<M> 
operator/(const int a, const ModAr<M> b) { return ModAr<M>(a)/b;}

template<int M> inline ModAr<M> 
operator==(const int a, const ModAr<M> b) { return b==a;}

template<int M> inline ModAr<M> 
operator!=(const int a, const ModAr<M> b) { return b!=a;}

// Define and set the static variables
template <int M>
int ModAr<M>::showmod = 0;   // show modulus when printing

template <int M>
ModAr<M>
ModAr<M>::operator^(const int e) const
{
   // Returns a^e.  Caution: the precedence of this operator does not follow
   // conventional mathematics.  Must use with parenthesis, as in a*(b^3).
   int ehere = e;
   int p=1;
   // checking if e==0 is faster than checking i<=e
   // if p==0, there is no need to go further
   for(; ehere&&p; ehere--)
	  p = (p*v) % M;  // mod at each step just to keep growth in check
   // (this is slower than it needs to be)
   return ModAr(p);
}

template <int M>
ModAr<M>
ModAr<M>::operator/(const ModAr &b)  const
{
   // mod operation is done in constructor
   if(b.v==0) {
	  exception zeroDivideException;
	  throw zeroDivideException;
	  //cout << "Division by zero not allowed\n";
	  //exit(-1);
   }
   int x,y,g;
   gcd(M,b.v,x,y,g);
   if(g != 1) {
	  exception singularNumberException;
	  throw singularNumberException;
	  //cerr << "Error: Attempting to divide by noninvertible number\n";
	  //exit(-1);
   }

   return ModAr(v*y);
}

template <int M>
ModAr<M>&
ModAr<M>:: operator/=(const ModAr &a) 
{
   if(a.v==0) {
	  exception zeroDivideException;
	  throw zeroDivideException;
	  // cout << "Division by zero not allowed\n";
	  // exit(-1);
   }
   int x,y,g;
   gcd(M,a.v,x,y,g);
   if(g != 1) {
	  exception singularNumberException;
	  throw singularNumberException;
	  // cerr << "Error: Attempting to divide by noninvertible number\n";
	  // exit(-1);
   }
   v = (v * y) % M;
   if(v < 0) v += M;
   return *this;
}


template <int M>
void ModAr<M>::gcd(const int m,const int b,int &x,int &y,int &g)
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

template <int M>
ostream& operator<<(ostream &os, const ModAr<M> &a) 
{
   if(ModAr<M>::getshowmod()) 
	  return os << a.getv() << '(' << M << ')';
   return os << a.getv();
}

#endif

/*
Local Variables:
compile-command: "g++ -c -g ModArnew.cc"
End:
*/

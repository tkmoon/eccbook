// GFNUM2m.h -- GF(2^m) numbers
// Created by Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef GFNUM2m_H
#define GFNUM2m_H
#include <iostream>
using namespace std;

#include "GF2.h"

enum outformat{vector,power};  // output options for GFNUM2m
class GFNUM2m;					// declare the class (for the externs below)

extern GFNUM2m ALPHA;			// set up a global alpha
extern GFNUM2m& A;				// and a reference to alpha for shorthand

class GFNUM2m {
protected:
   int v;				        // the value of the element (in vector form)

   static int gfm;				// number of elements in vector representation
   static int gfN;				// number of nonzero elments in field
   static outformat outtype;	// output type (vector or power representation)
   static int *v2p;				// convert from vector to power notation
								// (there is no exponential notation for zero)
   static int *p2v;				// convert from power to vector representation
   static void gferror(const char *s) {
	  cerr << s;
   } // print an error message
public:
   // static functions
   static void initgf(int m); // initialize multiplication array
   static void initgf(int m, unsigned int g); // initialize mult. array
   static void cleargf() {		// get rid of previous tables
	  delete [] p2v;
	  delete [] v2p;
	  p2v = 0;  v2p = 0;
   }
   static void setouttype(outformat ot) { outtype = ot;} // set the output type

   // Constructors
   GFNUM2m(const GFNUM2m& newnum) { v = newnum.v;}	// assignment 
   GFNUM2m(const int newv) { v = newv;};
                                 // pass in a number in vector representation
   GFNUM2m(void) { v = 0;};		 // default constructor 0 element)
   GFNUM2m(const GF2 gf) { v = gf.getv();};
   ~GFNUM2m() {};				// destructor

   // Functions to get member information
   int character() const { return 2; };
   int getv(void) const { return v; }; 
                                // return the vector representation
   int getp(void) const { return v2p[v]; };
                                // return the power representation
   static int getN(void) { return gfN; }; 
                                // return N, number of nonzero elts.
   static int getm(void) { return gfm; }; 
                                // return m (degree of vector sp.)
   static outformat getouttype() { return outtype; }; 
                                // return the output type
   
   // ***********************************************************************
   //                   the overloaded operators
   // ***********************************************************************

   GFNUM2m& operator+=(const GFNUM2m &a) { v ^= a.v;  return *this; }
   GFNUM2m& operator-=(const GFNUM2m &a) { v ^= a.v;  return *this; }
   GFNUM2m& operator*=(const GFNUM2m &a) {
	  if(a.v == 0 || v==0) v = 0;
	  else
		 v = p2v[(v2p[v]+a.v2p[a.v]) % gfN];
	  return *this;
   }
   GFNUM2m& operator/=(const GFNUM2m &a) {
	  if(a.v == 0) {
		 gferror("Division by zero in Galois field arithmetic\n");
	  }
	  else {
		 if(v == 0) return *this;
		 int x = v2p[v] - v2p[a.v];
		 if(x < 0) x += gfN;
		 v = p2v[x % gfN];
	  }
	  return *this;
   }
   GFNUM2m& operator^=(const int exp) {
	  if(v==0) {
		 if(exp == 0) {			//  make 0^0 = 1
			v = 1;
			return *this;
		 }
		 else {
			return *this;
		 }
	  }
	  if(exp < 0) {
		 v = p2v[gfN - ((v2p[v]*(-exp)) % gfN)];
	  }
	  else {
		 v = p2v[(v2p[v]*exp)% gfN];
	  }
	  return *this;
   }
	// unary -
   GFNUM2m operator-() const { return *this; }

   GFNUM2m operator+(const GFNUM2m &b) const {return GFNUM2m(v^b.v);}
   GFNUM2m operator-(const GFNUM2m &b) const {return GFNUM2m(v^b.v);}
   GFNUM2m operator*(const GFNUM2m &b) const {
	  GFNUM2m temp;
	  if(v == 0 || b.v == 0) return temp;
	  temp.v = p2v[(v2p[v] + v2p[b.v]) % gfN];
	  return temp;
   };
   GFNUM2m operator/(const GFNUM2m &b) const {
	  if(b.v == 0) {
		 gferror("Division by zero in Galois field arithmetic\n");
		 return 0;
	  }
	  else {
		 GFNUM2m temp;
		 if(v == 0) return temp;
		 int x = v2p[v] - v2p[b.v];
		 if(x < 0) x += gfN;
		 temp.v = p2v[x % gfN];
		 return temp;
	  }
   }

   // exponentiation operator: watch the precedence!!
   // const int e is changed to int e because there is no
   // way the variable location the user supplies for e is
   // going to be modified.  The function is more efficient
   // if it can modify its own copy of e.
   GFNUM2m operator^(const int e) const {
	  GFNUM2m temp;
	  if(v == 0) {
		 if(e == 0) {  // make 0^0 = 1
			temp.v = 1;
		 }
		 return temp;
	  }
	  if(v == 1) {
		 temp.v = 1;
		 return temp;
	  }
	  if(e < 0) {
		 temp.v = p2v[gfN- ((v2p[v]*(-e)) % gfN)];
	  }
	  else {
		 temp.v = p2v[v2p[v]*e % gfN];
	  }
	  return temp;
   }

   int operator==(const GFNUM2m &b) const {return (v == b.v);}
   int operator!=(const GFNUM2m &b) const {return (v != b.v);}

   friend ostream& operator<<(ostream &os, const GFNUM2m &a);
};

inline GFNUM2m operator+(const int a, const GFNUM2m b) { return b+a;}
inline GFNUM2m operator-(const int a, const GFNUM2m b) { return b+a;}
inline GFNUM2m operator*(const int a, const GFNUM2m b) { return b*a;}
inline GFNUM2m operator/(const int a, const GFNUM2m b) { return GFNUM2m(a)/b;}
inline GFNUM2m operator==(const int a, const GFNUM2m b) { return b==a;}
inline GFNUM2m operator!=(const int a, const GFNUM2m b) { return b!=a;}

#endif

/*
Local Variables:
compile-command: "g++ -c -g GFNUM2m.cc"
End:
*/

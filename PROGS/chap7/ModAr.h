// ******************************************************************
// MODAR.h -- MODulo ARithmetic
// Created by Todd K. Moon, with modifications by Stewart Weber
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef ModAr_H
#define ModAr_H

#include <iostream>
using namespace std;

class ModAr {
protected:
   int m, v;            // the modulo and value
   static int defaultm; // the default modulo value
   static int showmod;	// whether to display the modulo value
   void gcd(const int m,const int b,int &x,int &y,int &g) const;
           // Static member function used to compute inverses
public:
	// default constructor
   ModAr(int inv = 0, int inm= ModAr::defaultm) {
	  m = inm;
	  v = inv % m;  // Note: m should be >= 2: perhaps add test
	  if(v < 0) 
		 v += m;
   }
	// copy constructor
   ModAr(const ModAr& M) { v = M.v;  m = M.m;  }

	// destructor
   ~ModAr() {};	

   int getv() const { return v; };
   int getm() const { return m; };
   int character() const { return m; };
	// the int cast operator isn't included because it is identical
	// to the getv function
   int toint() { return (int)v; };  

   static void setdefaultm(int inm) { defaultm = inm; };
   static void setshowmod(int inm) { showmod = inm; };
   static int getdefaultm() { return defaultm; };
   
   // ***********************************************************************
   //                   the overloaded operators
   // ***********************************************************************

   ModAr& operator+=(const ModAr &a) {
	  v = (v + a.v) % m;
	  return *this;
   }
   ModAr& operator-=(const ModAr &a) {
	  v = (v - a.v) % m;
	  if(v < 0)  v += m;
	  return *this;
   }
   ModAr& operator*=(const ModAr &a) {
	  v = (v * a.v) % m;
	  return *this;
   }
   ModAr& operator/=(const ModAr &a);  // defined in ModAr.cc

	// unary -
   ModAr operator-() const { return ModAr(m - v,m); }

   ModAr operator+(const ModAr &b) const {return ModAr(v+b.v,m);}
   ModAr operator-(const ModAr &b) const {return ModAr(v-b.v,m);}
   ModAr operator*(const ModAr &b) const {return ModAr(v*b.v,m);}
   ModAr operator/(const ModAr &b) const; // defined in ModAr.cc

	// exponentiation operator: watch the precedence!!
	// const int e is changed to int e because there is no
	// way the variable location the user supplies for e is
	// going to be modified.  The function is more efficient
	// if it can modify its own copy of e.
   ModAr operator^(const int e) const; // defined in ModAr.cc

   int operator==(const ModAr &b) const {return (v == b.v) && (m == b.m);}
   int operator!=(const ModAr &b) const {return (v != b.v) || (m != b.m);}

   friend ostream& operator<<(ostream &os, const ModAr &a);
};

inline ModAr operator+(const int a, const ModAr b) { return b+a;}
inline ModAr operator-(const int a, const ModAr b) { return ModAr(a)-b;}
inline ModAr operator*(const int a, const ModAr b) { return b*a;}
inline ModAr operator/(const int a, const ModAr b) { return ModAr(a)/b;}
inline ModAr operator==(const int a, const ModAr b) { return b==a;}
inline ModAr operator!=(const int a, const ModAr b) { return b!=a;}

#endif

/*
Local Variables:
compile-command: "g++ -c -g ModAr.cc"
End:
*/

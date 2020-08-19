// GF2.h -- GF(2) numbers
// Created by Todd K. Moon

#ifndef GF2_H
#define GF2_H

#include <iostream>
using namespace std;
//#include <iostream>

class GF2 {
protected:
   unsigned char v;
   static void gferror(const char *s) {
	  cerr << s;
   } // print an error message
public:
   // Constructors
   GF2(const GF2& newnum) { v = newnum.v;}	// assignment 
   GF2(const int newv) { v = newv;};
                                 // pass in a number in vector representation
   GF2(void) { v = 0;};		 // default constructor 0 element)
   ~GF2() {};				// destructor

   // Functions to get member information
   int getv(void) const { return v; }; 
   // ***********************************************************************
   //                   the overloaded operators
   // ***********************************************************************

   GF2& operator+=(const GF2 &a) { v ^= a.v;  return *this; }
   GF2& operator-=(const GF2 &a) { v ^= a.v;  return *this; }
   GF2& operator*=(const GF2 &a) {
	  v &= a.v;
	  return *this;
   }
   GF2& operator/=(const GF2 &a) {
	  if(a.v == 0) {
		 gferror("Division by zero in Galois field arithmetic\n");
	  }
	  return *this;
   }
   GF2& operator^=(const int exp) {
	  return *this;
   }
	// unary -
   GF2 operator-() const { return *this; }

   GF2 operator+(const GF2 &b) const {return GF2(v^b.v);}
   GF2 operator-(const GF2 &b) const {return GF2(v^b.v);}
   GF2 operator*(const GF2 &b) const {
	  GF2 temp(v & b.v);
	  return temp;
   };
   GF2 operator/(const GF2 &b) const {
	  if(b.v == 0) {
		 gferror("Division by zero in Galois field arithmetic\n");
		 return 0;
	  }
	  else {
		 GF2 temp(v & b.v);
		 return temp;
	  }
   }

   // exponentiation operator: watch the precedence!!
   // const int e is changed to int e because there is no
   // way the variable location the user supplies for e is
   // going to be modified.  The function is more efficient
   // if it can modify its own copy of e.
   GF2 operator^(const int e) const {
	  return *this;
   }

   int operator==(const GF2 &b) const {return (v == b.v);}
   int operator!=(const GF2 &b) const {return (v != b.v);}

   friend ostream& operator<<(ostream &os, const GF2 &a) {
	  return os << int(a.v);
   }

};

inline GF2 operator+(const int a, const GF2 b) { return b+a;}
inline GF2 operator-(const int a, const GF2 b) { return b+a;}
inline GF2 operator*(const int a, const GF2 b) { return b*a;}
inline GF2 operator/(const int a, const GF2 b) { return GF2(a)/b;}
inline GF2 operator==(const int a, const GF2 b) { return b==a;}
inline GF2 operator!=(const int a, const GF2 b) { return b!=a;}

#endif

/*
Local Variables:
compile-command: "g++ -c -g GF2.cc"
End:
*/

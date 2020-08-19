// polynomialT.h -- templatized polynomial arithmetic
// Todd K. Moon, with modifications due to Stewart Weber

// This class uses a linked list to build a set of temporary
// polynomials, so that complicated expressions can be correctly parsed.

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

// changes since last version:  Added polynomialT<T>:: qualifier on lines
//  291, to make it work with new g++ version 4.2.2
#ifndef polynomialT_H
#define polynomialT_H

#include <cstdarg>
#include<iostream>
#include <assert.h>
using namespace std;

// The macro POLYC(TYPE,NAME,...) allows for easy 
// instatiation/initialization  of a polynomial with 
// a list of coefficients.  For example,
//
// POLYC(int,p1,{1,2,3,4});
//
// will create p1 as a polynomialT<int> type with the polynomial 
// 1 + 2x + 3x^2 + 4x^3
// in it.  The macro works by creating a polynomial,
// then creating an array initialized with the given list,
// then filling the polynomial with the list.
// However, the auxiliary array goes out of scope, so
// it does not occupy any futher space.
//
// This macro makes it convenient to initialize a polynomial array.
//
// Unfortunately, this macro does not work for the Microsoft compiler
// (as of this date) since their preprocessor does not appear
// to support macros with variable length arguements.
// For an example of how this can be replaced, see testpoly1.cc

#ifdef __GNUC__
#define POLYC(TYPE,NAME,vals...) polynomialT<TYPE> NAME; \
    { TYPE NAME##array[] = vals; \
      NAME.setc(NAME##array,sizeof(NAME##array)/sizeof(TYPE)-1); \
    }
// #define POLYC(TYPE,NAME,vals...) polynomialT<TYPE> NAME; \
//     { TYPE NAME##array[] = {vals}; \
//       NAME.setc(NAME##array,sizeof(NAME##array)/sizeof(TYPE)-1); \
//     }
#endif


template <class T>
class polynomialT {
//protected:
public:
   T *coeff;					// array of coefficients
   int degree;				  // degree (degree+1 elements are stored)
   void resize(int newdegree);
   char *varname;

   // the following static variables apply to all instances of this type
   static int decreasingprint; // Specifies order when printing: 
								//0=increasing (default), 1=decreasing
   static char defvarname[6];      // default variable name for printing

   static char beforeafterprint[3][6];
public:
   polynomialT(const T coeff0 = T(0), const char *vname=NULL); // default constructor
   polynomialT(int indegree, const T *incoeff, const char*vname=NULL);
								// constructor from an array of coeffs
   polynomialT(const polynomialT<T> &inpoly, const char *vname=NULL);
                                // copy constructor

   ~polynomialT();				// destructor

   void resizecopy(int newdegree);  // resize, copy old into new
   int getdegree(void) const { return degree; };
   polynomialT<T>& setc(int idx, const T &cval); // set a coefficient, 
   polynomialT<T>& setc(const T *cvals, int indegree);  
                                       // set indegree+1 coefficients 
   // setc always changes the degree as necessary, checking for 0 coefficients
   polynomialT<T>& buildspace(int indegree);
   // allocate space for a polynomial of at least this degree 
   // (even if some coefficients are zero)

   T* getarray(void ) { return coeff; }; // return array of coefficients


   // Overloaded operators
   // arithmetic {+, -, *, /, %, unary -}
   const polynomialT<T>& operator+(const polynomialT<T> &r) const;
   const polynomialT<T>& operator-(const polynomialT<T> &r) const;
   const polynomialT<T>& operator*(const polynomialT<T> &r) const;
   const polynomialT<T>& operator-(void) const;
   const polynomialT<T>& operator/(const polynomialT<T> &r) const;
   const polynomialT<T>& operator%(const polynomialT<T> &r) const;

   // assignment {=, +=, -=, *=, %= }
   const polynomialT<T>& operator=(const polynomialT<T> &r);
   const polynomialT<T>& operator+=(const polynomialT<T> &r);
   const polynomialT<T>& operator-=(const polynomialT<T> &r);
   const polynomialT<T>& operator*=(const polynomialT<T> &r);
   const polynomialT<T>& operator/=(const polynomialT<T> &r);
   const polynomialT<T>& operator%=(const polynomialT<T> &r);

   // division/remainder basic function
   const void divide(const polynomialT<T> &b) const;

   // Get results of last division or remainder
   static const polynomialT<T>& getlastquotient(void);
   static const polynomialT<T>& getlastremainder(void);

   // Shift operators (non cyclic shifts)
   const polynomialT<T>& operator>>(const int nshift) const; // divide by x^n
   const polynomialT<T>& operator<<(const int nshift) const ;//multiply by x^n
   const polynomialT<T>& operator>>=(const int nshift);	// divide by x^n
   const polynomialT<T>& operator<<=(const int nshift);	//multiply by x^n

   // indexing [] and evaluation ()
   T& operator[](int idx) const; // set a coefficient (in the range 0..deg)
   // i.e.,  poly[idx] = cval.
   // Be Careful: This function does not check the degree after assignment
   // However, it does ensure that idx satisfies 0 <= idx <= degree
   // Use function setc if you want to change the degree and check degree
   T& getcoeff(int i) const{assert(i >= 0 && i <= degree); return coeff[i]; }
   T operator()(const T &x) const;  // evaluate at x
   T evaluate(const T &x) const;    // evaluate at x
   
   // comparison {==, !=}
   int operator==(const polynomialT<T> &r) const;
   int operator!=(const polynomialT<T> &r) const;
   int operator==(const T &r) const;  // test for equality to a constant
   int operator!=(const T &r) const;  // test for inequality to a constant
   
   // A print member function so won't have to make << a friend,
   // but you may use    ostream << poly    for stream output
   ostream& printpolynomialT(ostream &os) const;
   ostream& printvarname(ostream &os) const;
   // prints the coefficient, possibly with a short string
   // before and after if the name matches that set in beforeafterprint
   ostream& printcoeff(ostream &os, const T &coeff) const;

   void setvarname(const char *newvarname);

   // set static variables governing all instances
   static void setprintdir(int dir); // 0=increasing, 1=decreasing
   static void setdefvarname(const char *newvarname);

   // set "before" and "after" strings used when the variable "name"
   // is printed
   static void setbeforeafterprint(const char *name,const char *before,
								   const char *after);


};

// output stream
template <class T> inline ostream&
operator<<(ostream &os, const polynomialT<T> &r) 
{
	r.printpolynomialT(os);
	return os;
}

// here are the operator special cases
template <class T> inline polynomialT<T>&
operator+(polynomialT<T> &l, const T &r) { return l + polynomialT<T>(r); }

template <class T> inline polynomialT<T>&
operator+(const T &l, polynomialT<T> &r) { return polynomialT<T>(l) + r; }

template <class T> inline polynomialT<T>&
operator-(const T &l, polynomialT<T> &r) { return polynomialT<T>(l) - r; }

template <class T> inline polynomialT<T>&
operator*(const T &l, polynomialT<T> &r) { return polynomialT<T>(l) * r; }

template <class T> inline const polynomialT<T>&
operator+=(polynomialT<T> &l, const T &r) { return l += polynomialT<T>(r); }

template <class T> inline const polynomialT<T>&
operator-=(polynomialT<T> &l, const T &r) { return l -= polynomialT<T>(r); }

template <class T> inline const polynomialT<T>&
operator*=(polynomialT<T> &l, const T &r) { return l *= polynomialT<T>(r); }


template <class T>
// class polytemp : public template <class T>  polynomialT<T>
// class polytemp :  template <class T> public polynomialT<T>
class polytemp : public polynomialT<T>
{
// protected:
public:
   // largest degree the polytemp can handle without
   // allocating more memory
   // remember that 'degree' is also a member of polytemp
   int maxdegree;
   
   // pointer to next available polytemp (makes linked list)
   polytemp<T>* next;
   
// The following static variables and functions are global to 
// all polytemps.

// here is a visual aid to show how the global variables are set up

//                    _________________________
//                   |                        
//  start_temp --->  | {polynomialT, maxdegree, _next_}______
//                   |_________________________      |
//                   |                             <- 
//                   | {polynomialT, maxdegree, _next_}______ 
//                   |_________________________      |
//                               ...                ... 
//                               ...                ...
//                    _________________________      |
//                   |                             <-     numpolytemp 
//  current_temp ->  | {polynomialT, maxdegree, _next_}______  (degree of list 
//                   |_________________________      |    up to here)
//                   |                             <- 
//  next_temp ---->  | {polynomialT, maxdegree, _next_}______ 
//                   |_________________________      |
//                               ...                ... 
//                               ...                ...
//                    _________________________      |
//  (last allocated  |                             <-     maxpolytemp
//   temppoly)       | {polynomialT, maxdegree, next=0}           (total degree
//                   |_________________________           of list)
// pointer to first polytemp in linked list
   static polytemp<T> *start_temp;
   
   // pointer to polytemp in current use, polytemps before
   // this one are being used, polytemps after this one
   // are available to be used
   static polytemp<T> *current_temp;
   
   // pointer to next available polytemp in linked list
   // this is equal to the 'next' member of current_temp
   static polytemp<T> *next_temp;
   
   // number of used polytemps (number of nodes in list from start_temp
   // to current_temp)
   static int numpolytemp;
   
   // total size of linked list (used and unused polytemps)
   static int maxpolytemp;
   
   // protected member functions
   void tempresize(int newdegree);
   
   
public:
   // constructor
   // the constructor for polynomialT doesn't care about degree
   // (it makes the degree 0 regardless) so a new polytemp
   // will need to have its size changed before use
   polytemp(void) : polynomialT<T>(T(0)), maxdegree(0), next(NULL) { };
   
   // static public member functions (available to anyone,
   // and global to all polytemps)
   
   // render all used polytemps as being unused
   static inline void resettemplist(void) {
	  next_temp = start_temp;
	  current_temp = NULL;
	  numpolytemp = 0;
   };

   // return pointer to an unused polytemp with a specific degree
   static polytemp<T>* gettemppoly(int degree);
   
   // print information about list size, etc.
   static void dumptemps(void);

   // get temporary quotient and remainder
   static polytemp<T>* gettempquot(int newdegree);
   static polytemp<T>* gettemprem(int newdegree);

//    static polytemp<T>* polytemp<T>:: gettempquot(int newdegree);
//    static polytemp<T>* polytemp<T>:: gettemprem(int newdegree);
// this version did not work on Microsoft V. 6

   // pointer to quotient
   static polytemp<T> *tempquotient;
   // pointer to remainder
   static polytemp<T> *temprem;
   void setdegree(int indegree) { if(indegree <= maxdegree) polynomialT<T>::degree = indegree;}

   const polytemp<T>& operator=(const polytemp<T> &r);
   const polytemp<T>& operator=(const polynomialT<T> &r);

};

// #define POLYC(TYPE,NAME,...) polynomialT<TYPE> NAME; \
//     { TYPE NAME##array[] = __VA_ARGS__; \
//       NAME.setc(NAME##array,sizeof(NAME##array)/sizeof(TYPE)-1); \
//     }

#endif

/*
Local Variables:
compile-command: "g++ -c -g polynomialT.cc"
End:
*/

// polynomialT.C -- templatized polynomial arithmetic
// Todd Moon, with modifications due to Stewart Weber

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include <assert.h>
#include "polynomialT.h"
extern "C" {
#include <string.h>  // for strlen
}

#ifndef MAX
#define MAX(a,b) (a>b?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) (a<b?(a):(b))
#endif

// changes since last version:  Added qualifier polynomialT<T>::
// on lines 926, 928, 929, 1056, 1064, 1066, to make it work with new g++.

template <class T> int polynomialT<T>::decreasingprint = 0;
template <class T> char polynomialT<T>::defvarname[6] = "x";
template <class T> char polynomialT<T>::beforeafterprint[3][6] = {{0},{0},{0}};

// default constructor
// same as constant assignment constructor
template <class T>
polynomialT<T>:: polynomialT(const T coeff0, const char *invarname) 
   : degree(0) 
{
   coeff = new T[1];
   coeff[0] = coeff0;
   
   varname = NULL;
   if(invarname) {
	  varname = new char[strlen(invarname)+1];
	  strcpy(varname,invarname);
   }
}


// normal constructor
template <class T>
polynomialT<T>:: polynomialT(int indegree, const T *incoeff, 
							 const char *invarname) 
{
   // check indegree
   for( ; (indegree>0) && (incoeff[indegree]==0); indegree--); 

   degree = indegree;
   coeff = new T[degree+1];
   for( ; indegree>=0; indegree--) // copy over coefficients
		coeff[indegree] = incoeff[indegree];
   polytemp<T>::resettemplist();

   varname = NULL;
   if(invarname) {
	  varname = new char[strlen(invarname)+1];
	  strcpy(varname,invarname);
   }
}


// copy constructor
template <class T>
polynomialT<T>:: polynomialT(const polynomialT<T> &inpoly, 
							 const char *invarname) 
{
   int i;

	i = degree = inpoly.degree;
	coeff = new T[degree+1];
	for( ; i>=0; i--)
		coeff[i] = inpoly.coeff[i];
	polytemp<T>::resettemplist();

	varname = NULL;
	if(inpoly.varname) {
	   varname = new char[strlen(inpoly.varname)+1];
	   strcpy(varname,inpoly.varname);
	}
}

template <class T>
polynomialT<T>::~polynomialT()
{ 
   delete [] coeff;
   if(varname) {
	  delete[] varname;
   }
}

// Protected member functions

// resize the polynomial (no copy)
template <class T> void
polynomialT<T>:: resize(const int newdegree)
{
   if(degree != newdegree) {
	  degree = newdegree;
	  delete [] coeff;
	  coeff = new T[newdegree+1];
   }
}

// resize and copy over the old 
// If the polynomial is made bigger, insert zeros in the 
// newly allocated spaces.  This makes things easier for 
// the functions that call resizecopy.

template <class T> void
polynomialT<T>:: resizecopy(const int newdegree)
{
   T *newcoeff;
   int i;
   T z(0);						// build a zero
   if(degree != newdegree) {
	  newcoeff = new T[newdegree+1];
	  for(i=newdegree; i>degree; i--)
		 newcoeff[i] = z;
	  for( ; i>=0; i--)
		 newcoeff[i] = coeff[i];
	  degree = newdegree;
	  delete [] coeff;
	  coeff = newcoeff;
   }
}


// Public member functions

// set a coefficient (can be used to change degree)
template <class T> polynomialT<T>&
polynomialT<T>:: setc(int idx, const T &cval) 
{
	int i;
	
	assert(idx>=0);
	
	if(idx <= degree) {
	   coeff[idx] = cval;
	   // fix poly if leading zeros
	   for(i=degree; (i>0) && (coeff[i]==0); i--);
	   if(i != degree) resizecopy(i);
	} else if(cval != 0) {
	   // note that idx>degree here
	   resizecopy(idx);
	   coeff[idx] = cval;
	}
	return *this;
}


// set a whole array of coefficients, setting degree based on the array
template <class T> polynomialT<T>&
polynomialT<T>:: setc(const T *cvals, int indegree) 
{
   int i;
   int d;
   assert(indegree>=0);
   // determine the actual degree
   for(d = indegree; (cvals[d] == 0) && (d > 0); d--) ;
   indegree = d;
   
   // make just enough room for the array
   resize(indegree);
   
   // copy over
   for(i = 0; i <= degree; i++) {
	  coeff[i] = cvals[i];
   }
   return *this;
}

// set a coefficient (can be used to change degree)
template <class T> polynomialT<T>&
polynomialT<T>:: buildspace(int indegree)
{
	int i;
	
	assert(indegree>=0);

	if(indegree > degree) {
	   resize(indegree);
	}
	return *this;
}


// Overloaded operators

// polynomial addition
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator+(const polynomialT<T> &r) const 
{
   polytemp<T> *out;
   int i;
   
   // There are 3 cases for deciding the size of
   // out (and 3 methods will be used to add), 
   // if l>r, then out is the size of 'l' 
   // (because we can be confident that 'l' doesn't
   // have any leading zeros).  If l<r, then out is
   // the size of 'r'.  If l==r, then we need to find
   // the highest nonzero coefficient of l+r to determine
	// the size.
   
   if(degree>r.degree) {
	  out = polytemp<T>::gettemppoly(degree);
	  for(i=degree; i>r.degree; i--)
		 out->coeff[i] = coeff[i];
	  for( ; i>=0; i--)
		 out->coeff[i] = coeff[i] + r.coeff[i];
   } else if(degree<r.degree) {
	  out = polytemp<T>::gettemppoly(r.degree);
	  for(i=r.degree; i>degree; i--)
		 out->coeff[i] = r.coeff[i];
	  for( ; i>=0; i--)
		 out->coeff[i] = coeff[i] + r.coeff[i];
   } else {
	  // degree==r.degree
	  // look for first nonzero coefficient
	  for(i=degree; (i>0) && (coeff[i]+r.coeff[i] == 0); i--);
	  out = polytemp<T>::gettemppoly(i);
	  for( ; i>=0; i--)
		 out->coeff[i] = coeff[i] + r.coeff[i];
   }
   
   return *(out);
} // overloaded +


// polynomial subtraction
// refer to operator+ for explanation of code
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator-(const polynomialT<T> &r) const 
{
   polytemp<T> *out;
   int i;

   if(degree>r.degree) {
	  out = polytemp<T>::gettemppoly(degree);
	  for(i=degree; i>r.degree; i--)
		 out->coeff[i] = coeff[i];
	  for( ; i>=0; i--)
		 out->coeff[i] = coeff[i] - r.coeff[i];
	  
   } else if(degree<r.degree) {
	  out = polytemp<T>::gettemppoly(r.degree);
	  
	  for(i=r.degree; i>degree; i--)
		 out->coeff[i] = -(r.coeff[i]);
	  for( ; i>=0; i--)
		 out->coeff[i] = coeff[i] - r.coeff[i];
	  
   } else {
	  // degree==r.degree
	  // look for first nonzero coefficient
	  for(i=degree; (i>0) && (coeff[i]-r.coeff[i] == 0); i--);
	  
	  out = polytemp<T>::gettemppoly(i);
	  
	  for( ; i>=0; i--)
		 out->coeff[i] = coeff[i] - r.coeff[i];
   }
   
   return *(out);
} // overloaded -


// polynomial multiplication
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator*(const polynomialT<T> &r) const
{
   polytemp<T> *out;
   T z(0);						// build a zero
   T endvalue = z;
   int i = degree + r.degree;
   int j, jdone;
   
   // Multiplication is just convolution.
   // Figure out degree of output by doing the
   // convolution of the endpoint first, then work
   // backwords.  The degree is the first nonzero coefficient
   
   // convolution without storage
   for( ; (i>=0) && (endvalue==0); i--) {
	  jdone = MIN(i, degree);
	  for(j=MAX(0, i-r.degree); j<=jdone; j++)
		 endvalue += (coeff[j]*r.coeff[i-j]);
   }
   out = polytemp<T>::gettemppoly(i+1);
   out->coeff[i+1] = endvalue;
   
   // finish convolution with storage
   for( ; i>=0; i--) {
	  out->coeff[i] = z;
	  
	  jdone = MIN(i, degree);
	  for(j=MAX(0, i-r.degree); j<=jdone; j++)
		 out->coeff[i] += (coeff[j]*r.coeff[i-j]);
   }
   return *out;
} // overloaded *

// polynomial division

// Since the operation of division produces both a quotient
// and a remainder, and sometimes both are wanted, both the quotient
// and remainder are saved in temporary variables.  They can then be
// retrieved using the functions getlastquotient and getlastremainder.

template <class T> const void 
polynomialT<T>:: divide(const polynomialT<T> &b) const
// compute this/b, with quotient and remainder
{
   int i,j,k;
   T z(0);						// build a zero
   if(degree < b.degree) {  // improper division
	  polytemp<T> *quot = polytemp<T>::gettempquot(0);
	  polytemp<T> *rem = polytemp<T>::gettemprem(degree);
	  quot->coeff[0] = z;
	  for(i = 0; i <= degree; i++) {
		 rem->coeff[i] = coeff[i];
	  }
   }
   else if(b.degree==0) {
	  T gninv = T(1)/b.coeff[0];
	  polytemp<T> *quot = polytemp<T>::gettempquot(degree);
	  polytemp<T> *rem = polytemp<T>::gettemprem(0);
	  for(i = 0; i <= degree; i++) {
		 (*quot)[i] = coeff[i] * gninv;
	  }
	  (*rem)[0] = z;
   }
   else {
	  int quotdeg = degree - b.degree;
	  int nr1 = b.degree-1;
	  T q, gninv;
	  int i,j,k;
	  polytemp<T> *quot = polytemp<T>::gettempquot(quotdeg);
	  polytemp<T> *rem = polytemp<T>::gettemprem(nr1);
	  // set up initial
	  gninv = T(1)/b.coeff[b.degree];
	  for(j = degree, i=nr1; i>=0; j--, i--) {
		 rem->coeff[i] = coeff[j];
	  }
	  for(k=0; k <= quotdeg; j--,k++) {  // work the rest of them
		 q = rem->coeff[nr1]*gninv;
		 quot->coeff[quotdeg-k] = q;
		 for(i = nr1; i>0; i--) {
			rem->coeff[i] = rem->coeff[i-1] - q*b.coeff[i];
		 }
		 rem->coeff[0] = coeff[j]-q*b.coeff[0];
	  }
	  // determine the actual degree of the remainder
	  rem->degree = 0;
	  for(i = nr1; i >=0; i--) {
		 if(rem->coeff[i] != 0) {
			rem->degree = i;
			break;
		 }
	  }
   }
}


template <class T> const polynomialT<T>& 
polynomialT<T>:: operator/(const polynomialT<T> &r) const
{
   divide(r);
   return *polytemp<T>::tempquotient;
}


template <class T> const polynomialT<T>& 
polynomialT<T>:: operator%(const polynomialT<T> &r) const
{  
   divide(r);
   return *polytemp<T>::temprem;
} // overloaded %


// polynomial div (/=)
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator/=(const polynomialT<T> &r) 
{
   (*this)=(*this)/r;
   polytemp<T>::resettemplist();
   return *this;
} // overloaded /=

// polynomial div (%=)
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator%=(const polynomialT<T> &r) 
{
   (*this)=(*this)%r;
   polytemp<T>::resettemplist();
   return *this;
} // overloaded %=

template <class T> const polynomialT<T>&
polynomialT<T>:: getlastquotient(void) {
   return *polytemp<T>::tempquotient;
}

template <class T> const polynomialT<T>&
polynomialT<T>:: getlastremainder(void) {
   return *polytemp<T>::temprem;
}


// unary -
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator-(void) const {
   polytemp<T>* out;
   int i;

   out = polytemp<T>::gettemppoly(degree);
   for(i=degree ; i>=0; i--)
	  out->coeff[i] = -(coeff[i]);
   return *(out);
} // overloaded unary -


// assignment
template <class T> const polynomialT<T>&
polynomialT<T>:: operator=(const polynomialT<T> &r) 
{
   int i = r.degree;

   if(this != &r) {
	  resize(r.degree);
	  
	  for( ; i>=0; i--)
		 coeff[i] = r.coeff[i];
	  
   }
   polytemp<T>::resettemplist();
   varname = NULL;
   if(r.varname) {
	  varname = new char[strlen(r.varname)+1];
	  strcpy(varname,r.varname);
   }
   return *this;
} // overloaded =


// polynomial increment (+=)
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator+=(const polynomialT<T> &r) 
{
   int i = r.degree;

	if(degree==r.degree) {
	   // look for first nonzero coefficient
	   for( ; (i>0) && (coeff[i]+r.coeff[i]==0); i--);
	   resizecopy(i);
	   
	} else if(degree<r.degree) resizecopy(r.degree);
	
	for( ; i>=0; i--) coeff[i] += r.coeff[i];
	
	polytemp<T>::resettemplist();
	return *this;
} // overloaded +=


// polynomial decrement (-=)
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator-=(const polynomialT<T> &r) 
{
   int i = r.degree;

   if(degree==r.degree) {
	  // look for first nonzero coefficient
	  for( ; (i>0) && (coeff[i]-r.coeff[i]==0); i--);
	  resizecopy(i);
	  
   } else if(degree<r.degree) resizecopy(r.degree);
   
   for( ; i>=0; i--) coeff[i] -= r.coeff[i];
   
   polytemp<T>::resettemplist();
   return *this;
} // overloaded -=


// polynomial scale (*=)
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator*=(const polynomialT<T> &r) 
{
   (*this)=(*this)*r;
   polytemp<T>::resettemplist();
   return *this;
} // overloaded *=


// divide by x -- reduce the degree (non-cyclic shift)
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator>>(const int nshift) const 
{
   polytemp<T> *out;
   T save;
   int i;
   T z(0);						// build a zero

   if(nshift == 0) {
	  return *this;
   }
   else if(nshift < 0) {		// shifting the other way
	  return (*this) << -nshift;
   }
   else {
	  if(nshift > degree) {
		 out = polytemp<T>::gettemppoly(0);
		 out->coeff[0] = z;
	  }
	  else {
		 out = polytemp<T>::gettemppoly(degree-nshift);
		 for(i=degree; i>nshift-1; i--) {
			out->coeff[i-nshift] = coeff[i];
		 }
	  }
   }
   return *out;
} // overloaded >>

// multiply by x -- increase the degree (non-cyclic shift)
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator<<(const int nshift) const 
{
   polytemp<T> *out;
   int i;
   T save;
   T z(0);						// build a zero

   if(nshift == 0) {			// no shift
	  return *this;
   }
   else if(nshift < 0) {		// shifting the other way
	  return (*this) >> -nshift;
   }
   else {
	  out = polytemp<T>::gettemppoly(degree+nshift);
	  for(i=degree; i>=0; i--) {
		 out->coeff[i+nshift] = coeff[i];
	  }
	  for(i = 0; i < nshift; i++) {
		 out->coeff[i] = z;
	  }
   }
   return *out;
} // overloaded <<

// multiply by x -- increase the degree (non-cyclic shift)
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator<<=(const int nshift) 
{
   int i;
   T save;
   T z(0);						// build a zero
   T *newcoeff;
   if(nshift == 0) {			// no shift
	  polytemp<T>::resettemplist();
	  return *this;
   }
   else if(nshift < 0) {		// shifting the other way
	  int degleft = degree - nshift;
	  if(degleft < 0) {			// nothing left
		 newcoeff = new T[1];
		 newcoeff[0] = z;
		 delete[] coeff;
		 coeff = newcoeff;
		 degree = 0;
	  }
	  else {					// stuff left
		 newcoeff = new T[degleft+1];
		 for(i = 0; i <= degleft; i++) {
			newcoeff[i] = coeff[i+nshift];
		 }
		 delete [] coeff;
		 coeff = newcoeff;
		 degree = degleft;
	  }
   }
   else {
	  newcoeff = new T[degree+nshift+1];
	  for(i=degree; i>=0; i--) {
		 newcoeff[i+nshift] = coeff[i];
	  }
	  for(i = 0; i < nshift; i++) {
		 newcoeff[i] = z;
	  }
	  degree = degree+nshift;
	  delete [] coeff;
	  coeff = newcoeff;
   }
   polytemp<T>::resettemplist();
   return *this;
} // overloaded <<

// divide by x -- decrease the degree (non-cyclic shift)
template <class T> const polynomialT<T>& 
polynomialT<T>:: operator>>=(const int nshift) 
{
   int i;
   T save;
   T z(0);						// build a zero
   T *newcoeff;
   if(nshift == 0) {			// no shift
	  polytemp<T>::resettemplist();
	  return *this;
   }
   else if(nshift < 0) {		// shifting the other way (increae the degree)
	  newcoeff = new T[degree+nshift+1];
	  for(i=degree; i>=0; i--) {
		 newcoeff[i+nshift] = coeff[i];
	  }
	  for(i = 0; i < nshift; i++) {
		 newcoeff[i] = z;
	  }
	  degree = degree+nshift;
	  delete [] coeff;
	  coeff = newcoeff;
   }
   else {
	  int degleft = degree - nshift;
	  if(degleft < 0) {			// nothing left
		 newcoeff = new T[1];
		 newcoeff[0] = z;
		 delete[] coeff;
		 coeff = newcoeff;
		 degree = 0;
	  }
	  else {					// stuff left
		 newcoeff = new T[degleft+1];
		 for(i = 0; i <= degleft; i++) {
			newcoeff[i] = coeff[i+nshift];
		 }
		 delete [] coeff;
		 coeff = newcoeff;
		 degree = degleft;
	  }
   }
   polytemp<T>::resettemplist();
   return *this;
} // overloaded <<


// indexing []
template <class T> T&
polynomialT<T>:: operator[](const int idx) const // const function
{
   assert(idx >= 0 && idx <= degree);
   return coeff[idx];
} // overloaded []


// evaluate the polynomial at the point 'x' (same as function evaluate)
template <class T> T
polynomialT<T>:: operator()(const T &x) const 
{
   T polyval = T(0);
   int i = degree;
   
   for(; i>=0; i--)
	  polyval = polyval*x + coeff[i];
   
   return polyval;
} // overloaded ()

// evaluate the polynomial at the point 'x'
template <class T> T
polynomialT<T>:: evaluate(const T &x) const 
{
   T polyval = T(0);
   int i = degree;
   
   for(; i>=0; i--)
	  polyval = polyval*x + coeff[i];
   
   return polyval;
} // evaluate

// compare (test if equal)
template <class T> int 
polynomialT<T>:: operator==(const polynomialT<T> &r) const
{
   int i = degree;
   
   if(degree != r.degree) return 0;
   
   for(; i>=0; i--)
	  if(coeff[i] != r.coeff[i]) return 0;
   
   return 1;
} // overloaded ==


// compare (test if not equal)
template <class T> int
polynomialT<T>:: operator!=(const polynomialT<T> &r) const 
{
	int i = degree;
	
	if(degree != r.degree) return 1;
	
	for(; i>=0; i--)
	   if(coeff[i] != r.coeff[i]) return 1;
	
	return 0;
} // overloaded !=


// test for equality to constant
template <class T> int 
polynomialT<T>:: operator==(const T &r) const 
{
   T z(0);						// build a zero
   // if any higher terms are nonzero, cannot be equal
   for(int i = 1; i <= degree; i++) {
	  if(coeff[i] != z) return 0;
   }
   if(coeff[0] != r) return 0;
   return 1;
} // overloaded ==


// test for inequality to constant
template <class T> int 
polynomialT<T>:: operator!=(const T &r) const 
{
   T z(0);						// build a zero

   // if any higher terms are nonzero, then it is not equal
   for(int i = 1; i <= degree; i++) {
	  if(coeff[i] != z) return 1;
   }
   if(coeff[0] == r) return 0;
   return 1;
} // overloaded ==



template <class T> ostream&
polynomialT<T>:: printvarname(ostream &os) const 
{
   char *vn;
   if(varname) {
	  vn = varname;
   }
   else {
	  vn = polynomialT<T>::defvarname;
   }
   os << vn;
   return os;
}

template <class T> ostream&
polynomialT<T>:: printcoeff(ostream &os, const T &coeff) const 
{
   char *vn;
   if(varname) {
	  vn = varname;
   }
   else {
	  vn = polynomialT<T>::defvarname;
   }
   if(!strcmp(vn,beforeafterprint[0]))
	  os << beforeafterprint[1];
   os << coeff;
   if(!strcmp(vn,beforeafterprint[0]))
	  os << beforeafterprint[2];
   return os;
}


// print function
template <class T> ostream&
polynomialT<T>:: printpolynomialT(ostream &os) const 
{
   int i;

   if(degree==0) {
	  os << coeff[0];
	  return os;
   }

   if(polynomialT<T>::decreasingprint == 0) { // print in increasing order
	  // do first two terms to take care of the possiblity
	  // that they are the 0 and 1 degree coefficients
	  
	  // find start (first nonzero term)
	  for(i=0; coeff[i]==0 && i < degree; i++);
	  
	  if((i==0) || ((i>0) && (coeff[i]!=1))) printcoeff(os,coeff[i]);
	  
	  if(i>0) printvarname(os);
	  if(i>1) os << "^" << i;
	  
	  // is there only 1 term?
	  if(i>=degree) return os;
	  
	  // find next non-zero term
	  for(i++; coeff[i]==0 && i < degree; i++);
	  
	  os << " + ";
	  if(coeff[i]!=1) printcoeff(os,coeff[i]);
	  printvarname(os);
	  if(i>1) os << "^" << i;
	  if(i>=degree) return os;
	  
	  // output remaining polynomial
	  for(i++; i<degree; i++) {
		 if(coeff[i] == 1) {
			os << " + " ;
			printvarname(os);
			os << "^" << i;
		 }
		 else if(coeff[i] != 0) {
			os << " + ";  printcoeff(os,coeff[i]);
			printvarname(os);
			os << "^" << i;
		 }
	  }
	  // output highest degree term always (in case it is 0)
	  os << " + ";
	  if(coeff[i] != 1) {
		 printcoeff(os,coeff[i]);
		 printvarname(os);
		 os << "^" << degree;
	  }
	  else {
		 printvarname(os);
		 os << "^" << degree;
	  }
   }
   else { // print in decreasing order
	  // always print the first term (in case it is 0)
	  if(coeff[degree] != 1) printcoeff(os,coeff[degree]);
	  printvarname(os);
	  if(degree > 1) os << "^" << degree;
	  
	  // find next
	  for(i = degree-1; i >= 0; i--) {
		 if(coeff[i] != 0) {
			os << " + ";
			if(coeff[i] != 1 || i==0) printcoeff(os,coeff[i]);
			if(i > 0) {
			   printvarname(os);
			}
			if(i > 1) os << "^" << i;
		 }
	  }
   }
   return os;
}

template <class T> void 
polynomialT<T>::setprintdir(const int dir)
{ decreasingprint=dir;
}


template <class T> void
polynomialT<T>::setdefvarname(const char *newvarname)
{ if(strlen(newvarname) < 5)
	 strcpy(polynomialT<T>::defvarname,newvarname);
}

template <class T> void
polynomialT<T>::setbeforeafterprint(const char *name,const char *before,
									const char *after)
{
   if(strlen(name) < 6) 
	 strcpy(polynomialT<T>::beforeafterprint[0],name);
   if(strlen(before) < 6)
	 strcpy(polynomialT<T>::beforeafterprint[1],before);
   if(strlen(after) < 6)
	 strcpy(polynomialT<T>::beforeafterprint[2],after);
}


template <class T> void
polynomialT<T>::setvarname(const char *newvarname)
{ 
   if(varname) {
	  delete[] varname;
   }
   varname = new char[strlen(newvarname)+1];
   strcpy(varname,newvarname);
}

// polytemp class functions

// Initialize the protected static variables
template <class T> polytemp<T>*
polytemp<T>:: start_temp = NULL;

template <class T> polytemp<T>*
polytemp<T>:: current_temp = NULL;

template <class T> polytemp<T>*
polytemp<T>:: next_temp = NULL;

template <class T> int 
polytemp<T>:: numpolytemp = 0;

template <class T> int 
polytemp<T>:: maxpolytemp = 0;

// Initialize variables for quotient and remainder
template <class T> polytemp<T>* polytemp<T>:: tempquotient = NULL;
template <class T> polytemp<T>* polytemp<T>:: temprem = NULL;

// polytemp class protected member functions

// resize coeff array and set new maxsize
template <class T> void 
polytemp<T>:: tempresize(const int newdegree) {
   polynomialT<T>::degree = newdegree;
   maxdegree = newdegree;
   delete [] polynomialT<T>::coeff;
   polynomialT<T>::coeff = new T[newdegree+1];
}


// polytemp class public static functions

// return pointer to an unused polytemp with a specific size
template <class T> polytemp<T>*
polytemp<T>:: gettemppoly(const int degree) 
{
   polytemp<T>* temp;

   // We need to do one of four things:
   // 1) Return pointer to next available polytemp. 
   // 2) Resize next available polytemp, then return it.
   // 3) Create a new poly temp and add it on to the list.
   // 4) Create the first node in the list.
   
   // This will be the order that these cases are addressed
   // because on average, the higher number cases are less
   // likely (case 4 is only addressed once).

   // case 1 and 2
   if(polytemp<T>::numpolytemp < polytemp<T>::maxpolytemp) {
	   
	  polytemp<T>::numpolytemp++;
	  
	  temp = polytemp<T>::current_temp = polytemp<T>::next_temp;
	  polytemp<T>::next_temp = temp->next;
	  
	  temp->degree = degree;
	  
	  // case 2 (case 1 if this test is false)
	  if(degree > temp->maxdegree) temp->tempresize(degree);
	  
	  return temp;
   }
   
   // case 3 and 4
   // (note: For case 3 and 4, we know that next_temp = NULL
   // coming into this function.  Since we are adding onto
   // the list at this point, next_temp is going to remain
   // NULL.)
   
   polytemp<T>::numpolytemp++;
   polytemp<T>::maxpolytemp++;
   
   temp = new polytemp<T>;
   temp->tempresize(degree);
   
   // case 3 
   // (note: maxpolytemp was just incremented, so
   // to test if it used to be > 0, test if it is > 1)
   if(polytemp<T>::maxpolytemp > 1) 
	  polytemp<T>::current_temp->next = temp;
   
   // case 4 (start list for first time)
   else
	  polytemp<T>::start_temp = temp;
   
   polytemp<T>::current_temp = temp;
   
   return temp;
}


// print information about list size, etc.
template <class T> void
polytemp<T>:: dumptemps() 
{

   polytemp<T>* temp = polytemp<T>::start_temp;
   int counttemps = 0;
   
   cout<<"Temporary information:"<<endl
	   <<"maxpolytemp="<<polytemp<T>::maxpolytemp<<endl;
   
   for( ; temp!=NULL; temp=temp->next, counttemps++)
	  cout<<"ctr="<<counttemps// <<" address: "<<temp
		  <<" maxdegree="<<temp->maxdegree<<": "
		  <<*temp<<endl;
   
   if(counttemps != polytemp<T>::maxpolytemp)
	  cout<<"Warning: queue length does not match number allocated"
		  <<endl;
   if(polytemp<T>::tempquotient) { // have a quotient
	  cout << "Temp quotient: " << *polytemp<T>::tempquotient << endl;
   }
   if(polytemp<T>::temprem) { // have a quotient
	  cout << "Temp remainder: " << *polytemp<T>::temprem << endl;
   }
}

template <class T> polytemp<T>*
polytemp<T>:: gettempquot(int newdegree) 
{
   if(polytemp<T>::tempquotient == NULL) {
	  polytemp<T>::tempquotient = new polytemp<T>;
   }
   if(newdegree > polytemp<T>::tempquotient->maxdegree)
	  polytemp<T>::tempquotient->tempresize(newdegree);
   polytemp<T>::tempquotient->degree = newdegree;
   return polytemp<T>::tempquotient;
}

template <class T> polytemp<T>*
polytemp<T>:: gettemprem(int newdegree) 
{
   if(polytemp<T>::temprem == NULL) {
	  polytemp<T>::temprem = new polytemp<T>;
   }
   if(newdegree > polytemp<T>::temprem->maxdegree)
	  polytemp<T>::temprem->tempresize(newdegree);
   polytemp<T>::temprem->degree = newdegree;
   return polytemp<T>::temprem;
}


// assignment
template <class T> const polytemp<T>&
polytemp<T>:: operator=(const polytemp<T> &r) 
{
   if(this != &r) {
	  polynomialT<T>::degree = r.degree;
	  for(int i = r.degree; i>=0; i--)
		 polynomialT<T>::coeff[i] = r.coeff[i];
   }
   // polytemp<T>::resettemplist(); // don't do this: other's might grab this
   return *this;
} // overloaded =

template <class T> const polytemp<T>&
polytemp<T>:: operator=(const polynomialT<T> &r) 
{
   if(this != &r) {
	  polynomialT<T>::degree = r.getdegree();
	  for(int i = r.getdegree(); i>=0; i--)
		 polynomialT<T>::coeff[i] = r[i];
   }
   // polytemp<T>::resettemplist(); // don't do this: other's might grab this
   return *this;
} // overloaded =


/*
Local Variables:
compile-command: "g++ -c -g polynomialT.cc"
End:
*/

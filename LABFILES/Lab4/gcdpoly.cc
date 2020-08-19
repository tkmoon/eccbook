// gcdpoly.cc --- Implement a gcd algorithm for a general
// polynomial type.  Also, implement a gcd for polynomials with real
// coefficients, truncating coefficients as necessary to avoid
// roundoff problems.
//
// Copyright 2019 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include <math.h>
#include "polynomialT.h"

// declare the gcd function, setting up a default value for the rdeg
template <class T> void
gcd(const polynomialT<T> &a, const polynomialT<T> &b, 
	polynomialT<T> &g,
	polynomialT<T> &s, polynomialT<T> &t, int rdeg=0);

// Function definition: gcd
template <class T>  void
gcd(const polynomialT<T> &a, const polynomialT<T> &b, polynomialT<T> &g,
	polynomialT<T> &s, polynomialT<T> &t, int rdeg)
{
   polynomialT<T> qi;  // quotient
   polynomialT<T> ri, rim1, rim2, si, sim1, sim2, ti, tim1, tim2;
   polynomialT<T> *rip, *rim1p, *rim2p, *sip, *sim1p, *sim2p;
   polynomialT<T> *tip, *tim1p, *tim2p;
   polynomialT<T> *temp;
   T norm;

   // Fill in the blanks...
}


// static void chop(polynomialT<double> *p, double eps);

// // This is a specializtion for doubles, since it has to handle
// // the roundoff more carefully
// template <> void
// gcd(const polynomialT<double> &a, const polynomialT<double> &b, 
// 	polynomialT<double> &g,	polynomialT<double> &s, polynomialT<double> &t,
// 	int rdeg)
// {
//    polynomialT<double> qi;  // quotient
//    polynomialT<double> ri, rim1, rim2, si, sim1, sim2, ti, tim1, tim2;
//    polynomialT<double> *rip, *rim1p, *rim2p, *sip, *sim1p, *sim2p;
//    polynomialT<double> *tip, *tim1p, *tim2p;
//    polynomialT<double> *temp;
//    double norm;

//    rim2 = a;
//    rim1 = b;
//    sim2 = 1;  sim1 = 0;
//    tim2 = 0;  tim1 = 1;

//    rip = &ri;  rim1p = &rim1;  rim2p = &rim2;
//    sip = &si;  sim1p = &sim1;  sim2p = &sim2;
//    tip = &ti;  tim1p = &tim1;  tim2p = &tim2;

//    do {
// 	  qi = *rim2p / *rim1p;
// 	  *rip = *rim2p - *rim1p*qi;
// 	  chop(rip,1e-10);
// 	  *sip = *sim2p - *sim1p*qi;
// 	  chop(sip,1e-10);
// 	  *tip = *tim2p - *tim1p*qi;
// 	  chop(tip,1e-10);
// 	  temp = rim2p; rim2p = rim1p;  rim1p = rip;  rip = temp;
// 	  temp = sim2p; sim2p = sim1p;  sim1p = sip;  sip = temp;
// 	  temp = tim2p; tim2p = tim1p;  tim1p = tip;  tip = temp;
//    }
//    while(rim1p->getdegree() != 0);
//    s = *sim2p;   t = *tim2p;   g = *rim2p;
//    norm = 1/g[g.getdegree()];
//    s *= norm;
//    g *= norm;
//    t *= norm;
// }

// void chop(polynomialT<double> *p, double eps)
// {
//    int i;
//    int newdegree=0;
//    int done = 0;
//    for(i = p->getdegree(); i>= 0; i--) {
// 	  if(fabs((*p)[i]) < eps) {
// 		 (*p)[i] = 0;
// 		 if(i== p->getdegree()) {
// 			newdegree = i;  // set that new degree is necessary
// 		 }
// 	  }
//    }
//    if(newdegree) {
// 	  for(i = newdegree-1; i > 0; i--) {
// 		 if((*p)[i] != 0) {
// 			p->resizecopy(i);
// 			done = 1;
// 			break;
// 		 }
// 	  }
// 	  if(!done) p->resizecopy(0);
//    }
// }

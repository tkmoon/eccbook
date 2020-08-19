//
//
//  Program: bchdesigner.cc --- find binary primitive BCH codes

//  Todd K. Moon
//  Utah State University
//
//  Date:  March 31, 2004

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only


#include <iostream>
using namespace std;
#include "GFNUM2m.h"
#include "polynomialT.cc"

extern "C" {
#include <stdlib.h>
}

// template class polynomialT<GFNUM2m>;

int main(int argc, char **argv)
{

   if(argc==1) {
	  cerr << "\nUsage: " << argv[0] << "[-n n] [-t t] [-b b] [-M] [-p] \n\n";
	  cerr << "Prints BCH designs for primitive binary BCH codes\n\n";
	  cerr << "-n -- code length.  Must be 2^m-1. default:15\n";
	  cerr << "-t -- design correction capability.  default:2\n";
      cerr << "-b -- specify starting exponent.  default:1 (narrow sense)\n";
  	  cerr << "-M -- search over all b to find one with largest\n";
	  cerr << "      dimension and longest consecutive sequence of roots\n";
	  cerr << "-p -- print the minimal polynomials which are the factors\n";
	  return(-1);
   }
   
   int n = 15;					// code length
   int t = 2;					// correction capability
   int best=0;					// find the best b
   int b = 1;
   int printM = 0;				// print minimal polynomials

   int i;

   for(i = 1; i < argc; i++) {
	  if(!strcmp(argv[i],"-n"))
		 n = atoi(argv[++i]);
	  if(!strcmp(argv[i],"-t"))
		 t = atoi(argv[++i]);
	  if(!strcmp(argv[i],"-b"))
		 b = atoi(argv[++i]);
	  if(!strcmp(argv[i],"-M"))
		 best = 1;
	  if(!strcmp(argv[i],"-p"))
		 printM = 1;
   }
   

   int k;						// code dimension
   int q = 2;					// binary codes.
   
   int m=0;						// used to define field
   for(i = 0; i < 32; i++) {
	  if((1<<i) == n+1) {
		 m = i;
		 break;
	  }
   }
   if(m == 0) {
	  cerr << "Error: cannot find field for this n.\n";
	  exit(-1);
   }
   GFNUM2m::initgf(m);			// initialize the field

   int delta = 2*t+1;
   int startb, endb;
   if(best == 0) {
	  startb = b;
	  endb = b;
   }
   else {
	  startb = 0;
	  endb = n-1;
   }

   int a;

   POLYC(GFNUM2m,l1,{0,1});		// linear term x + 0
   polynomialT<GFNUM2m> p;		// minimal polynomial
   polynomialT<GFNUM2m> g(1);	// generator polynomial

   int *expvals = new int[n];
   int ex;						// exponent value (modulo n)
   int bestrunlength = 0;
   int bestk = 0;
   polynomialT<GFNUM2m> gbest;
   int bbest,rlbest;

   for(b = startb; b <= endb; b++) {
	  g = GFNUM2m(1);
	  // cout << "b=" << b << endl;
	  for(i = 0; i < n+1; i++) expvals[i] = 0;
	  // roots: A^b,A^(b+1), ... A^(b+delta-2): total of 2t consecutive values
	  for(i = b; i <= b+delta-2; i++) {
		 ex = i % n;				//  wrap around
		 // cout << "i=" << i << "  ex=" << ex << endl;
		 if(expvals[ex]) {		// this one already; done
			continue;
		 }
		 // take this one and all its conjugates
		 expvals[ex] = 1;
		 a = ex;
		 l1[0] = A^ex;
		 p = l1;					// first factor of minimal polynomial
		 while(1) {
			a = (a*q) % n;
			// cout << "a=" << a << " ";
			if(a != ex) {
			   expvals[a] = 1;
			   l1[0] = A^a;
			   p = p * l1;
			}
			else {
			   break;
			}
		 }
		 // cout << endl;
		 if(printM) {
			cout << "Minimal polynomial=" << p << endl;
		 }
		 g *= p;					// factor in generator
	  }
	  // degree of generator = n-k
	  k = n - g.getdegree();
	  // find the longest run of consecutive roots
	  int runlength = 0, maxrunlength = 0;
	  for(i = 0; i < n; i++) {
		 if(expvals[i] == 1) {
			runlength++;
		 }
		 else {
			if(runlength) {		// a run was already started
			   if(runlength > maxrunlength) {
				  maxrunlength = runlength;
			   }
			   runlength = 0;
			}
		 }
	  }
	  if(runlength) {	  // if still in a run, see if it wraps around
		 for(i = 0; i < n; i++) {
			if(expvals[i] == 1) {
			   runlength++;
			}
			else {
			   if(runlength) {		// a run was already started
				  if(runlength > maxrunlength) {
					 maxrunlength = runlength;
				  }
				  runlength = 0;
			   }
			   break;				// this time, stop at the end of a run
			}
		 }
	  }
	  if(maxrunlength > bestrunlength) {
		 bestrunlength = maxrunlength;
	  }
	  if(k > bestk) {
		 bestk = k;
		 gbest = g;
		 bbest = b;
		 rlbest = maxrunlength;
	  }
   }  // end loop over b

   cout << "n=" << n << " t(des)=" << t << " b=" << bbest << 
	  " k=" << bestk << "  t(actual)="  << 	  (double)rlbest/2.0 << endl;
   cout << "g=" << gbest << endl;

}


// this stuff is stuck here where hopefully it won't be too readily 
// discovered, since it provides a solution to a lab exercise

// Define static variables in GFNUM2m
int *GFNUM2m::p2v = 0;			// convert exponent to vector
int *GFNUM2m::v2p = 0;			// a list of elements to convert from
								// vector to exponential notation
int GFNUM2m::gfm = 0;			// vector size of field element
int GFNUM2m::gfN = 0;			// number of nonzero-elements in field
outformat GFNUM2m::outtype = GFpower;
								// default to exponential output
// ALPHA elements from the field
GFNUM2m ALPHA;                  // define the element alpha
GFNUM2m& A= ALPHA;              // and a shorthand reference to it


void GFNUM2m::initgf(int m)
{
   //do the initialization by using only the size of the field
   // m in GF(2^m).
   // A fixed set of primitive polynomials is used.
   // Octal:
   unsigned int g[] = {1,1,7,013, 023, 045, 0103, 0211, 0435, 01021, 02011,
					   04005, 010123, 020033, 042103, 0100003};
   if(m>sizeof(g)/sizeof(unsigned int)) {
	  cerr << "Error: must specify connection polynomial for m" << endl;
   }
   GFNUM2m::initgf(m,g[m]);
}

// initgf: (1) Build the v2p and p2v tables
//         (2) set the global variable ALPHA
//         (3) set the static member variables gfm and gfN
void GFNUM2m::initgf(int m,unsigned int g)
// m -- GF(2^m)
// g -- generator polynomial, bits represent the coefficients:
// e.g. g = 0x13 = 1 0011 = D^4 + D + 1
{
   int i,j;

   if(m > sizeof(unsigned int)*8) {	// too many bits!
	  cerr << "Error: Degree too large in GFNUM2m" << endl;
   }
	  
   ALPHA.v = 2;					// set up alpha element
   gfm = m;
   gfN = (1<<m)-1;				// gfN = number of nonzero field elements,
								// gfN = 2^n -1

   if(v2p) delete[] v2p;		// delete any prior stuff
   if(p2v) delete[] p2v;
   
   v2p = new int[gfN+1];  		// table to convert vector to power form
   p2v = new int[gfN+1];		// tabel to convert power to vector form

   // convert g to an integer
   int gint=0;
   int mask1 = 1<<(m-1);
   int bitout;
   gint = g;
   v2p[0] = gfN;				// no conversion of this element

   int lfsr = 1;				// initial load of lfsr representation: alpha^0=1
   for(i = 1; i <= gfN; i++) {
	  v2p[lfsr] = i-1;
	  bitout = ((lfsr&mask1)!=0);
	  lfsr = ((lfsr<<1) ^ bitout*gint)&gfN;
   }

   p2v[0] = 1;
   for(i = 1;i <= gfN; i++) {
	  p2v[v2p[i]] = i;
   }
   p2v[gfN] = 1;
}


ostream& operator<<(ostream& s,const GFNUM2m& arg)
{
   if(arg.getouttype() == GFpower) {
	  if(arg.v < 2) return s << arg.v;
	  else {
		 int e = arg.v2p[arg.v];
		 if(e == 1) return s << "A";
		 else return s << "A^" << e;
	  }
   }
   else {						// vector (numeric) output
	  return s << arg.v;
   }
}


  
/*
Local Variables:
compile-command: "g++ -o bchdesigner -g bchdesigner.cc"
End:
*/



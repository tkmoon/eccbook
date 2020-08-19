//  Program: testGS.cc --- test the GS decoder and multivar poly arithmetic
//  This program actually goes through a decoding process for a RS code

//  Todd K. Moon, Feb 7, 2004


// Select one of the following two lines
#define GFTYPE
//#define MODARTYPE


#ifdef MODARTYPE
#include "ModAr.h"
#define TYPE ModAr
#endif
#ifdef GFTYPE
#include "GFNUM2m.h"
#define TYPE GFNUM2m
#endif


#include <math.h>

using namespace std;

#include "polynomialT.cc"
#include "rothruck.h"

// create instantiations of the polynomial class 
template class polynomialT<TYPE>;
template class polytemp<TYPE>;
template class polynomialT<polynomialT<TYPE> >;


int computetm(int n,int k,int m);
int computeLm(int n,int k, int m);

polynomialT<polynomialT<GFNUM2m> >
kotter(int n,int L,GFNUM2m *xi, GFNUM2m *yi,int *mi, int *wdeg,int rlex,
	   int m1);
GFNUM2m evaluate(polynomialT<polynomialT<GFNUM2m> > &Q,GFNUM2m a,
				 GFNUM2m b);
TYPE computeD(int r,int s,const polynomialT<polynomialT<TYPE> > &Q,
			  TYPE a, TYPE b);
int hammdist(const polynomialT<TYPE> &p, const polynomialT<TYPE> &r, 
			 const TYPE *ss, int n);

int main()
{
#ifdef MODARTYPE
   TYPE::setdefaultm(5);		// work over GF19
   TYPE::setshowmod(0);		// don't show modulus
#endif
#ifdef GFTYPE
   GFNUM2m::initgf(4,0x13);  // 1 0011 = 1+d+d^4
#endif

   int n = 15;
   int k = 7;
   int d = n-k+1;
   int t = int(floor((d-1.)/2.));
   int m1 = 2;					// interpolation multiplicity
   
   cout << "n=" << n << "  k=" << k << "  d=" << d << "  t0=" << t <<
	  "  m=" << m1 << "  tm=" << computetm(n,k,m1) << 
	  "  Lm=" << computeLm(n,k,m1) << endl;

   POLYC(GFNUM2m,m,{A,A^2,A^3,A^4,A^5,A^6,A^7});
   cout << "message polynomial=" << m << endl;
   GFNUM2m cc[n];				// code vector
   GFNUM2m sc[n];				// support set
   int supset1 = 0, supset2 = 14; // support set of code : 1,...A^14
   int i,j;
   for(j=0, i = supset1; i <= supset2; i++, j++) {
	  cc[j] = m(A^i);
	  sc[j] = A^i;
   }
   polynomialT<GFNUM2m> c(n-1,cc);
   polynomialT<GFNUM2m> ss(n-1,sc);
   cout << "code polynomial=" << c << endl;
   // build error 
   int nerr = 4;
   int errloc[n];
   GFNUM2m errval[n];
   errloc[0] = 1;
   errloc[1] = 3;
   errloc[2] = 5;
   errloc[3] = 7;
   errloc[4] = 9;
   // errloc[5] = 11;
   errval[0] = A^2;
   errval[1] = A^3;
   errval[2] = A^4;
   errval[3] = A^5;
   errval[4] = A^6;
   // errval[5] = A^7;
   GFNUM2m ec[n];
   for(i = 0; i < nerr; i++) {
	  ec[ errloc[i]] = errval[i];
   }
   polynomialT<GFNUM2m> e(n-1,ec);
   cout << "e=" << e << endl;
   polynomialT<GFNUM2m> r(c+e);	
   cout << "received polynomial=" << r << endl;
   int wdeg[2] = {1,k-1};		// weighted degree
   int rlex = 1;				// rlex=1 for rlex; rlex=0 for lex

   polynomialT<polynomialT<GFNUM2m> > Q; // this will be the interpolating poly
   Q.setvarname("y");			// set the way Q is printed
   Q.setbeforeafterprint("y","(",")");

   // if mi=NULL, then the same multiplicity is used by every point
   int *mi = NULL;

   int Lm = computeLm(n,k,m1);

   Q = kotter(n,Lm,ss.getarray(),r.getarray(),mi, wdeg,rlex,m1);
   cout << "interpolating polynomial=" << Q << endl;
   // check the results
   for(i = 0; i < n; i++) {
	  int C = (m1+1)*m1/2;
	  for(int rs = 0; rs < C; rs++) {
		 int s1 = rs % m1;
		 int r1 = rs / m1;
		 cout <<"Q("<<ss[i]<<"," << r[i] << ")_(" << r1 << "," << s1<< ")=" << 
			computeD(r1,s1,Q,ss[i],r[i]) << endl;
	  }
   }
   
//   cout << "Q[0]=" << Q[0] << endl;
//   cout << "Q[1]=" << Q[1] << endl;

//    polynomialT<GFNUM2m> P1(Q[1]);
//    polynomialT<GFNUM2m> P0(Q[0]);
//    polynomialT<GFNUM2m> f(P0 / P1);
//    polynomialT<GFNUM2m> rn(P0 % P1);
//    cout << "P1=" << P1 << "  P0 =" << P0 << "  f=" << f << "  r=" << rn << endl;
   rpolynode *rptr,*rptr1;		// pointer to list of polynomials
   rptr = rothruck(Q,k-1);		// find the y-roots of Q

   cout << "y-roots of Q(x,y): " << endl;
   for(rptr1 = rptr; rptr1 != NULL; rptr1 = rptr1->next) {
	  cout << rptr1->f << "  Hamming dist=" << hammdist(rptr1->f,r,sc,n) <<
		 endl;
   }
}


int hammdist(const polynomialT<TYPE> &p, const polynomialT<TYPE> &r, 
			 const TYPE *ss, int n)
{
   TYPE codesym;
   int hdist = 0;
   int j;
   for(j=0; j <= r.getdegree(); j++) {
	  codesym = p(ss[j]);
	  if(codesym != r[j]) hdist++;
   }
   for(  ; j < n; j++) {
	  codesym = p(ss[j]);
	  if(codesym != 0) hdist++;
   }
   return hdist;
}
	  
   

#include "kotter.cc"  
#include "rothruck.cc"
#include "computetm.cc"
#include "computeLm.cc"

/*
Local Variables:
compile-command: "g++ -o testGS3 -Wno-deprecated -g testGS3.cc ModAr.cc polynomialT.cc GFNUM2m.cc"
End:
*/




#include <math.h>
#include <iostream>
using namespace std;

#include "bittype.h"


double uran(void)
{
   return rand()/(double)RAND_MAX;
}

static int graniset = 0;
static double grangset;
double gran(void)
{
   double rsq, v1,
	  v2, fac;

   if(!graniset) {
	  graniset = 1;
	  do {
		 v1 = 2*(rand()/(double)RAND_MAX) - 1;
		 v2 = 2*(rand()/(double)RAND_MAX) - 1;
		 rsq = v1*v1 + v2*v2;
	  } while(rsq > 1 || rsq == 0);
	  fac = sqrt(-2*log(rsq)/rsq);
	  grangset = v1*fac;
	  return v2*fac;
   }
   else {
	  graniset = 0;
	  return grangset;
   }
}


//-----------------------------------------------------------------------

void bpskmodbits(BITTYPE *x, double *s, double Ecsqrt, int N)
// take the bits in x and BPSK modulate into s
//  0 --> -sqrt(Ec)     1 --> sqrt(Ec)
{
   for(int i = 0; i < N; i++) {
	  s[i] = Ecsqrt*(2*x[i] - 1);
   }
}

void randnoise(double *n, double sigma, int N)
// make AWGN noise vector 
{
   for(int i = 0; i < N; i++) {
	  n[i] = sigma*gran();
   }
}

void addnoise(double *s, double *n, double *y, int N)
// y = s + n
{
   for(int i = 0; i < N; i++) {
	  y[i] = s[i] + n[i];
   }
}

void addscalenoise(double *s, double *n, double sigma, double *y, int N)
// y = s + sigma*n
{
   for(int i = 0; i < N; i++) {
	  y[i] = s[i] + sigma*n[i];
   }
}


void randbits(BITTYPE *u, int N)
{
   for(int i = 0; i < N; i++) {
	  if(uran() > 0.5) u[i] = 1;
	  else u[i] = 0;
   }
}


/*
Local Variables:
compile-command: "clang++ -c randstuff.cc"
End:
*/

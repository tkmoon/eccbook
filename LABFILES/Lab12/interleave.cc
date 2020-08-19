// interleave.cc -- A random interleaver
// Todd K. Moon

// This is meant to be a very elementary interleaver, 
// operating simply on the basis of random selection.
// Interleavers requirless less memory (and less setup) 
// could be made (e.g., using a 2-d array).


#include "interleave.h"
#include <iostream>
using namespace std;
extern "C" {
#include <stdlib.h>
}
#include <math.h>
double uran(void);

interleave::interleave(int in_size, unsigned int seed)
{
   int i,j,k,r,numleft;
   double u;
   size = in_size;
   pi = new int[size];
   piinv = new int[size];
   double rm;
   if(seed) srand(seed); else srand(1);

   // generate a random permutation.  See Knuth, v. 2, p. 145
   for(i = 0; i < size; i++) {
   	  pi[i] = i;
   };
   for(j = 0; j < size; j++) {
   	  u = uran();
   	  k = floor(j*u);
   	  pi[j] = pi[k];
   	  pi[k] = j;
   	  piinv[pi[k]] = k;
   }
   for(j = 0; j < size; j++) {
	  piinv[pi[j]] = j;
   }
}

void
interleave::Pi(const double *in, double *out)
{
   int i;
   for(i = 0; i < size; i++) {
	  out[i] = in[pi[i]];
   }
}

void
interleave::Pi(const double **in, double **out, int nrow)
{
   int i,j;
   for(i = 0; i < size; i++) {
	  for(j = 0; j < nrow; j++) {
		 out[i][j] = in[pi[i]][j];
	  }
   }
}

void
interleave::Pi(const unsigned char *in, unsigned char *out)
{
   int i;
   for(i = 0; i < size; i++) {
	  out[i] = in[pi[i]];
   }
}

void
interleave::Piinv(const double *in, double *out)
{
   int i;
   for(i = 0; i < size; i++) {
	  out[i] = in[piinv[i]];
   }
}

void
interleave::Piinv(const double **in, double **out, int nrow)
{
   int i,j;
   for(i = 0; i < size; i++) {
	  for(j = 0; j < nrow; j++) {
		 out[i][j] = in[piinv[i]][j];
	  }
   }
}

void
interleave::PiinvTimesoverlay(const double **in, double **out, int nrow)
{
   int i,j;
   for(i = 0; i < size; i++) {
	  for(j = 0; j < nrow; j++) {
		 out[i][j] *= in[piinv[i]][j];
	  }
   }
}

void
interleave::Piinv(const unsigned char *in, unsigned char *out)
{
   int i;
   for(i = 0; i < size; i++) {
	  out[i] = in[piinv[i]];
   }
}


/*
Local Variables:
compile-command: "g++ -c interleave.cc"
End:
*/

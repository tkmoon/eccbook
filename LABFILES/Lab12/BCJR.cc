// BCJR.cc -- Compute the MAP bit decisions using the BCJR algorithm
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "BCJR.h"
#include <math.h>
#include <iostream>
using namespace std;
#include "BPSKmodvec.h"


BCJR::BCJR(BinConv &encoder, int in_L, double insigma2)
{
   L = in_L;
   sigma2 = insigma2;			// must be set by function call
   Q = (1<<encoder.nu);
   n = encoder.n;
   k = encoder.k;
   numbranches = (1<<k);
   // allocate space for alpha, beta, and the posterior probabilities
   CALLOCMATRIX(alpha,double,(L+1),Q);
   CALLOCMATRIX(beta,double,(L+1),Q);
   CALLOCMATRIX(Pp,double,L,numbranches);

   buildprev(encoder);
   buildoutputmat(encoder);
}



double **
BCJR::MAP1(const double **r, const double **prior, int finalstate)
// Compute the MAP probability for u, using the entire r vector at
// each step and using the entire gamma probability (not just the
// extrinsic part)
{
   int i,k,p,q;
   unsigned int in;
   int inp;

   alphabeta(r,prior,finalstate); // compute alpha and beta
   // Find the log likelihood ratios
   // Transition: (p,q):  S[1][p] = q for transition due to input 1
   // Transition: (p,q):  S[0][p] = q for transition due to input 0
   double sum,sum2;
   for(k = 0; k < L; k++) {
	  sum2 = 0;
	  for(inp = 0; inp < numbranches; inp++) {
		 sum = 0;
		 for(i = 0; i < Q; i++) {
			p = i;  q = S[inp][i];
			sum += alpha[k][p]*gammakpq(k,p,q,inp,prior[k][inp],r[k])*
			   beta[k+1][q];
		 }
		 Pp[k][inp] = sum;
		 sum2 += sum;
	  }
	  // normalize across all the inputs
	  for(inp = 0; inp < numbranches; inp++) {
		 Pp[k][inp] /= sum2;
	  }
   }

   return Pp;
}

double **
BCJR::MAP2(const double **r,const double **prior, int finalstate)
// Compute the MAP probability for u, using only the nonsystematic
// part of r.   This is used for computing the extrinsic probabilities

{
   int i,k,p,q;
   unsigned int in,inp;

   // Fill in the blanks

   return Pp;
}

void
BCJR::alphabeta(const double **r,const double **prior,int finalstate)
// Compute alpha and beta using the forward-backward algorithm.
// Set finalstate<0 if uniformly distributed.
// Otherwise, set finalstate to the desired value
{
   int i,p,q;
   double normalpha, normbeta, gamma;
   unsigned int in;
   double Ppq;					// transition probability
   // initialize alpha[0][*] 
   for(i = 0; i < Q; i++) {
	  alpha[0][i] = 0;
   }
   alpha[0][0] = 1;
   // initialize beta[0][*]
   if(finalstate>=0) { // initialize for one final state
	  for(i = 0; i < Q; i++) {
		 beta[L][i] = 0;
	  }
	  beta[L][finalstate] = 1;
   }
   else { // initialize all uniformly
	  for(i = 0; i < Q; i++) {
		 beta[L][i] = 1./Q;
	  }
   }

   // Forward pass: compute the alphas
   // Fill in the blanks ...
	  // normalize the alpha
      // Fillin the blanks ...
   // Backward pass: compute the betas

   // Fill in the blanks  ...
	  // normalize the beta
      // Fill in the blanks ...

  
}

double
BCJR::gammakpq(int k, int p, int q, int inp, double Ppq, const double *r)
// Compute the transition probability; assumes that (p,q) is a 
// valid state transition
// 
{
   double *a;					// output selection
   int j;

   a = outputmat[p][inp];
   double t, sum = 0;
   for(j = 0; j < n; j++) {
	  t = r[j] - a[j];
	  sum += t*t;
   }
   return Ppq*exp(-sum/(2*sigma2));
}

double
BCJR::gammakpq2(int k, int p, int q, int inp, const double *r)
// Compute the transition probability; assumes that (p,q) is a valid
// state transition
// 
// This version uses only the "nonsystematic" part of r
// and is used for computing the extrinsic probabilities
{
   double *a;					// output selection
   int j;

   a = outputmat[p][inp];
   double t, sum = 0;
   for(j = 1; j < n; j++) {
	  t = r[j] - a[j];
	  sum += t*t;
   }
   return exp(-sum/(2*sigma2));
}

		 
void 
BCJR::buildprev(BinConv &encoder)
{
   unsigned int savestate = encoder.getstate();

   CALLOCMATRIX(prevstate,unsigned int, Q, numbranches);
   CALLOCMATRIX(inputfrom,unsigned int, Q, numbranches);
   CALLOCMATRIX(S,unsigned int,numbranches,Q);
   
   // first build the nextstate 
   unsigned int **nextstate;
   CALLOCMATRIX(nextstate,unsigned int,Q,numbranches);
   unsigned char* ins = new unsigned char[k];
   unsigned int state;
   unsigned int inp;
   int i;
   unsigned int nextst;
   int *nfrom = new int[Q];

   for(state = 0; state < Q; state++) {
	  for(inp = 0; inp < numbranches; inp++)  {
		 encoder.setstate(state);
		 // convert inp to array
		 for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
		 encoder.encode(ins);
		 nextst = encoder.getstate();
		 prevstate[nextst][nfrom[nextst]] = state;
		 inputfrom[nextst][nfrom[nextst]] = inp;
		 nfrom[nextst]++;
		 S[inp][state] = nextst;
	  }
   }
	  
   for(state = 0; state < Q; state++) { // for each state
	  for(inp = 0; inp < numbranches; inp++)  {
		 encoder.setstate(state);
		 // convert inp to array
		 for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
		 encoder.encode(ins);
		 nextstate[state][inp] = encoder.getstate();
	  }
   }
   for(state = 0; state < Q; state++) {
	  unsigned int *outs;
	  cout << "state=" << state << ": ";
	  for(inp=0; inp < numbranches; inp++) {
		 cout << nextstate[state][inp] << " ";
		 encoder.setstate(state);
		 for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
		 outs = encoder.encode(ins);
		 cout << "(";
		 for(i = 0; i < n; i++) {
			cout << int(outs[i]);
		 }
		 cout << ") ";
	  }
	  cout << endl;
   }
   encoder.setstate(savestate);


   // print the prevstate table
//    cout << "fromstates: " << endl;
//    for(state = 0; state < Q; state++) {
// 	  cout << "state=" << state << ": ";
// 	  for(inp=0; inp < (1<<k); inp++) {
// 		 cout << prevstate[state][inp] << " (";
// 		 cout << inputfrom[state][inp] << ") ";
// 	  }
// 	  cout << endl;
//    }
   delete[] nfrom;
   delete [] ins;
   encoder.setstate(savestate);
   FREEMATRIX(nextstate);
}


void 
BCJR::buildoutputmat(BinConv & encoder)
// builds the lookup from [state][input] to output array,
{
   unsigned int savestate = encoder.getstate();

   // outputmat[state][input][outputnum]
   CALLOCTENSOR(outputmat,double,Q,numbranches,n);

   unsigned char* ins = new unsigned char[k];
   unsigned int state;
   unsigned int inp;
   int i,j;
   unsigned int *out;
   unsigned int outint;
   double *modout;
   BPSKmodvec modulator(n);

   for(state = 0; state < Q; state++) { // for each state
	  for(inp = 0; inp < numbranches; inp++)  {
		 encoder.setstate(state);
		 // convert inp to array
		 for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
		 out = encoder.encode(ins);
		 modout = modulator.mod(out);
		 for(j = 0; j < n; j++) { 
			outputmat[state][inp][j] = modout[j];
		 }
	  }
   }
   cout << "outputs: " << endl;
   for(state = 0; state < Q; state++) {
	  cout << "state=" << state << ": ";
	  for(inp=0; inp < numbranches; inp++) {
		 for(j = 0; j < n; j++) {
			cout << int(outputmat[state][inp][j]);
		 }
		 cout << " ";
	  }
	  cout << endl;
   }
   encoder.setstate(savestate);
   delete [] ins;
}		 

/*
Local Variables:
compile-command: "g++ -c -g BCJR.cc"
End:
*/

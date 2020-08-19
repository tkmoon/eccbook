// BCJR.h -- Compute the MAP bit decisions using the BCJR algorithm
// Todd K. Moon
// Feb 17, 2003

#ifndef BCJR_H
#define BCJR_H
#include "BinConv.h"
#include "matalloc.h"

class BCJR {
protected:

   int nu;						// constraint length
   int Q;						// Q = 2^nu; number of states
   int L;						// total number of trellis stages
   int n;						// dimension of received signal vector
   double sigma2;				// noise variance in each component
   int N1trans;					// number of state transitions where x=1
   unsigned int **S;			// (p,q): S[1][p]=q or S[0][p] = q
   unsigned int **prevstate;    // [Q][Q] list of previous states
   unsigned int **inputfrom;	// [Q][numbranches] -- inputs
   double ***outputmat;			// [Q][numbranches][n] -- outputs
   int k;						// number of inputs
   int numbranches;				// number of branches = 2^numinputs

   void buildprev(BinConv &encoder); // build information about trellis
   void buildoutputmat(BinConv & encoder);

   double gammakpq(int k, int p, int q, int inp, double Ppq, const double *r);
   // Compute the transition probability;
   // Assumes that (p,q) is a valid
   // state transition.
   // inp is the input at state p, 
   // Ppq is the probability of the input that causes the transition p to q
   double gammakpq2(int k, int p, int q, int i,const double *r);
   // This version is used for computing the "extrinsic probabilities."
   // It uses only the "nonsystematic" part of r
   // and does not use the prior probabilities.
public:
   // provide access to the forward/backward data (for testing purposes)
   double **alpha;				// [L+1][Q] forward probabilities
   double **beta;				// [L+1][Q] backward probabilities
   double **Pp;					// [L][numbranch] posterior bit probability


   BCJR(BinConv  &encoder, int in_L, double insigma2=0);
   ~BCJR() {FREEMATRIX(S); FREEMATRIX(inputfrom); FREETENSOR(outputmat,Q);
   FREEMATRIX(alpha);  FREEMATRIX(beta);  delete[] Pp;
   };

   void setsigma2(double in_sigma2) { sigma2 = in_sigma2; }
   // set the noise variance for the observations

   void alphabeta(const double **r,const double **prior, int finalstate=0);
   // compute the forward and backward passes
   // r[i][n] -- received sequence
   // prior[numbranch] -- the prior probabilities at time k for each input
   // e.g. prior[0] = P(u_k=0), prior[1] = P(u_k=1)
   double **MAP1(const double **r, const double **prior, int finalstate=0);
   // Compute the LLR using the entire r vector (not systematic)
   // r[i][n] -- received sequence
   // set finalstate<0 for 
   // uniform distribution, or set finalstate to actual state
   // returns a pointer to the posterior probability
   double **MAP2(const double **r, const double **prior, int finalstate=0);
   // Compute the LLR using nonsystematic part and systematic parts separately
   // Returns a pointer to the extrinsic part of the probability
   // which resides in Pp
};

#endif

/*
Local Variables:
compile-command: "g++ -c -g BCJR.cc"
End:
*/

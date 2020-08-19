// ldpcdec.h --- LDPC decoder class declarations
// Todd K. Moon

#ifndef LDPCDEC_H
#define LDPCDEC_H
#include <string>
#include "matalloc.h"
#include <iostream>
using namespace std;


class LDPCDEC {
public:  // for debugging purposes
   int decodetype;				// 1=probability decoder, 2=log likelihood

   // variables for reprsenting the matrix:
   int N;						// block length
   int K;						// message length
   int M;						// N-K -- redundancy
   int *Nmlen;					// [M] lengths of row weight vectors
   int **Nm;					// [M][rowwt]
                                //  set of bits that participate in check m
   int *Mnlen;					// [N] lengths of column weight vectors
   int **Mn;					// [N][colwt]
                                // set of checks in which bit n participates
   int maxcolwt;				// maximum weight of columns
   int maxrowwt;				// maximum weight of rows

   unsigned char *c;			// decoded bit values

   // variables for the probability decoder:
   int *na;						// [N] # of elements above this in column
   double *pn;					// [N] channel posterior probabilities
   double *deltaq;				// [maxrowwt] temporary row info holding deltaq
   double **r1;					// [maxcolwt][N]  r1[m][n]
   double **r0;					// [maxcolwt][N]  r0[m][n]
   double **q0;					// [maxcolwt][N]  q0[m][n]
   double **q1;					// [maxcolwt][N]  q1[m][n]
   double **deltar;				// [maxcolwt][N]  deltar[m][n]
   double *q0p;					// [N] -- pseudopriors
   double *q1p;					// [N] -- pseudopriors

   // variables for the log likelihood decoder
   double *rn;					// soft inputs
   double *lambda;				// log likelihoods
   double **u;					// check messages
   double **oldu;



   void allocdecodedat(void);	// allocate the decoder memory
   void freedecodedat(void);	// free the decoder memory

   // channel parameters

   double sigma2;				// channel variance
   double a;					// coded signal amplitude = sqrt(E_c)

   void probdecodeinit(double *y);	// initialize the prob. decoder
   void logdecodeinit(double *y);	// initialize the log like decodoer
   void probdecode();				// probability decoder
   void lldecode();					// log like decoder
public:
   int printstuff;				// set for printing things out
   LDPCDEC(const string fname, int offset=0, int decodetype=1); 
   // constructor --- read from file
   //  Some definitions use base index=1, while some uses base index=0.
   // Use offset=1 when reading a file with base index=1
   // Use offset=0 when reading a file with base index=0
   // decodtype: 0=probability, 1=log likelihood
   ~LDPCDEC() {
	  freedecodedat(); 
	  delete[] Mnlen;
	  delete[] Nmlen;
	  FREEMATRIX(Mn);
	  FREEMATRIX(Nm);
   };
   int decode(double *pn, int printstuff, int maxnumloop, int & numloops);
   void printsparse(const string& name, double **sparsemat);
   void printsparseA(void);
   void setsigma2(double ins2) { sigma2 = ins2; }
   void setsigamp(double ina) { a = ina; }

   void meastoprobAWGN(double *y);
   // convert measured data to a probability.  In this case, for an AWGN
   // with variance sigma2 and channel modulation amplitude a
};

#endif
/*
Local Variables:
compile-command: "g++ -c -g ldpcdec.cc"
End:
*/

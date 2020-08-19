// ldpcdecoder.h --- LDPC decoder class declarations
// Todd K. Moon

#ifndef LDPCDECODER_H
#define LDPCDECODER_H
#include "matalloc.h"
#include <string>
#include <fstream>

#define DOLPSTUFF
#ifdef DOLPSTUFF  // include code for linear programming decoding
extern "C" {
#include <stdio.h>
#include <glpk.h>
}
#include <vector>
#endif

using namespace std;

#ifndef MAX
#define MAX(a,b) (a) >= (b) ? (a) : (b)
#endif

// Define bit masks for different kinds of decoders
#define DOPROBDECODE 1			// 0
#define DOLOGLIKEDECODE 2		// 1
#define DOMINSUMDECODE 4		// 2
#define DOPHIDECODE 8			// 3
#define DOQPHIDECODE 16			// 4
#define DOMINSUMCORRECTDECODE 32 // 5
#define DOAMINSTARDECODE 64		 // 6
#define DORCBPDECODE 128		 // 7
#define DOBITFLIP1 256			 // 8
#define GALABITFLIP 512			 // 9
#define WEIGHTEDBITFLIP 1024	 // 10
#define MODWEIGHTEDBITFLIP 2048	 // 11
#define GDBITFLIP 4096			 // 12
#define MULTIGDBITFLIP 8192		 // 13
#define DC1 16384				 // 14
#define DMBP 32768				 // 15
#define LPDECODE 65536           // 16

extern std::string decodernameslist[];

class LDPCDECODER {
public:
   unsigned long int decodetype;  // bit mask of decoders defined above
   int numdecodetype;
   int *decodetypelist;		  // list of requested decoder types
   int *decodetypenum;
   //   std::string *decodenames;		  // points to name of type of decoder

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

   
   // A pool of data used by all decoders
   unsigned char **c;			// decoded bits (one for each decoder type)
   double *problikedata;		// the input probability or likelihood
   unsigned int *na;			// [N] # of elements above this in column
   double *rowdata;				// from row, for for/back processing
   double *forwardprod;			// table for for/back processing
   double *backwardprod;		// table for for/back processing
   double *forbackdataout;		// output for for/back processing
   double *outproblikedata1;	// output probability or likelihood
   double *outproblikedata0;	// output probability or likelihood (prob decod)
   double **mtondata1;			// [N][maxcolwt] decoder m to n
   double **mtondata0;			// [N][maxcolwt] decoder m to n (prob decod)
   double **ntomdata1;			// [maxcolwt][N]
   double **ntomdata0;			// [maxcolwt][N] (prob decode)

   unsigned long int *numloops;	// number of loops used by decoder
   // numloops = number of decoding iterations actually used
   // if numloops is a positive integer on input, then
   // the decoder does exactly that many iterations.  
   // if numloops is 0 on input, then the decoder
   // proceeds until parities all check, or until maxnumloop iterations

   // variables for the probability decoder:
   double *pn;					// [N] channel posterior probabilities
   double *deltaq;				// [maxrowwt] temporary row info holding deltaq
   double **Pmton1;					// [N][maxcolwt] Pmton1[m][n]
   double **Pmton0;					// [N][maxcolwt]  Pmton0[n][m]
   double **q0;					// [maxcolwt][N]  q0[m][n]
   double **Pntom0;					// [maxcolwt][N]  Pntom1[m][n]
   double **Pntom1;					// [maxcolwt][N]  Pntom1[m][n]
   double **deltar;				// [maxcolwt][N]  deltar[m][n]
   double *Pnout0;					// [N] -- pseudopriors
   double *Pnout1;					// [N] -- pseudopriors
   double *poutprods;
   double *poutprods0, *poutprods1; // leave-one-out product results [maxcolwt]
   double *forwardprod0, *backwardprod0, *forwardprod1, *backwardprod1;


   // variables for the log likelihood decoder
   double *Lc;					// Channel log likelihood
   double *tanhstuff;			// tanh information for this column
   double *poutprobsLL;           // leave-one-out sum probs
   double *Lcout;				// output posterior log
   double **Lntom;  // messages variable to check
   double **Lmton;  // messages check to variables
   double *poutsum;
   double *forwardsum0, *backwardsum0;

   // Variables for min sum decoder
   double *MSLc;				// Channel log likelihood
   double *MSLcout;
   char **alphantom;			// signs:  use 0 to indicate +, 1 to indicate -
   // then can compute using parity
   double **betantom;			// magnitudes
   double **MSLmton;			// var to check message
   double *betarow;
   char *alpharow;
   double *poutmin;
   double *forwardmin0;
   double *backwardmin0;

   // Variables for the phi decoder
   double *phiLc;				// channel log likelihood
   double *phiLcout;			// output likelihood
   double **phiLmton;

   // variables for bit flipping
   unsigned char *synd;
   unsigned int *flist;

   
   void allocdecodedat(void);	// allocate the decoder memory
   void freedecodedat(void);	// free the decoder memory

   // channel parameters

   double sigma2;				// channel variance
   double a;					// coded signal amplitude
   // if a is negative, then a bit 0 is represented by a negative value
   //  abs(a) = sqrt(E_c)
   int a_sign;					// sign of a.

   
   int printstuff;				// set for printing things out
   LDPCDECODER(std::string fname, int offset=0, unsigned long int decodetype=1); 
   // constructor --- read from file
   //  Some definitions use base index=1, while some uses base index=0.
   // Use offset=1 when reading a file with base index=1
   // Use offset=0 when reading a file with base index=0
   // decodtype: 0=probability, 1=log likelihood
   ~LDPCDECODER() { freedecodedat(); 
               delete[] Mnlen;
			   delete[] Nmlen;
			   FREEMATRIX(Mn);
			   FREEMATRIX(Nm);
			   delete [] decodetypelist;
			   delete [] decodetypenum;
			   };
   unsigned long int decode(double *y, int printstuff,
							unsigned long int maxnumloop,
							unsigned long int *numloops, double *durations);
   void printsparse(std::string name, double **sparsemat);
   void printsparseMp1(std::string name, double **sparsemat);
   void printsparsechar(std::string name, unsigned char **sparsemat);
   void printsparsechar(std::string name, char **sparsemat);
   void printsparsetranspose(std::string name, double **sparsemat);
   void printsparseA(void);
   void setsigma2(double ins2) { sigma2 = ins2; }
   void setsigamp(double ina) { a = ina;
	  if(a >= 0) a_sign = 1; else a_sign = -1;
}
   int checkparity(unsigned char *);  // check parity. c = (alleged) codeword

   int decodetypetonum(unsigned long int decodetype);

   void meastoprobAWGN(double *y);
   // convert measured data to a probability.  In this case, for an AWGN
   // with variance sigma2 and channel modulation amplitude a

   void forbackprobs(double *p, int m,double *poutprobs);
   // take an array of data and do the leave-one-out probabilities,
   // return result in poutprobs, and the overall product in pall
   void forbackprobs2(double b0, double *p0, double b1, double *p1, int m,
								double *poutprods0,double &pall0,
					  double *poutprods1,double &pall1);
   // same as forbackprobs, but works on 2 arrays of data, and
   // includes an overall factor to capture multiplication by pn[n] and 1-pn[n]
   void forbacksum(double b0, double *p0, int m,
				   double *poutsum0,double &psum0);
   // leave-one-out sum
   void forbackmin(double *p0, int m,double *poutmin0);
   // leave-one-out min

   // variables and functions for the quantized phi decoder
   int Mquantphi;				// number of quantiziation levels
   double xquantphimax;			// maximum value of x
   double phiquantdelta;			// quantization interval
   double *qphitable;			// table of quantized values
   double *phibetarow;			// row of phi(beta) data
   double *phioutphisum;		// phi sum from for/back
   double phi(double x);		// compute the phi function
   double phiq(double x);		// compute the quantized phi function
   double boxplus(double x, double y);  // compute the box plus function
   double rcbpboxplus(double x, double y);  // compute the box plus function
   void buildqphitable(void);	// build the table of quantized phi function
   void phidecodeinit(double *y);
   int phidecode(unsigned char *c, unsigned long int maxnumloop,
				 unsigned long int &numloops);
   int qphidecode(unsigned char *c, unsigned long int maxnumloop,
				  unsigned long int &numloops);

   void probdecodeinit(double *y);	// initialize the prob. decoder
   int probdecode(unsigned char *c, unsigned long int maxnumloop,
				  unsigned long int &numloops);

   void lldecodeinit(double *y);	// initialize the log like decodoer
   int lldecode(unsigned char *c, unsigned long int maxnumloop,
				  unsigned long int &numloops);

   void minsumdecodeinit(double *y);// initialize min sum decoder
   int minsumdecode(unsigned char *c, unsigned long int maxnumloop,
				  unsigned long int &numloops);


   // variables and functions for the min sum with correction decoder
   double *MScorrectLc;		// Channel log likelihood
   double *MScorrectLcout;
   double **MScorrectLmton;	// var to check message
   double scorrect(double x,double y);
   double stildecorrect(double x,double y);
   void minsumcorrectdecodeinit(double *y);
   int minsumcorrectdecode(unsigned char *c,unsigned long int maxnumloop,
							  unsigned long int &numloops);
   void forbackmincorrect(double *p0, int m,double *poutmin0);

   // variables and functions for the a-min* decoder
   double *aminstarLc;		// Channel log likelihood
   double *aminstarLcout;
   double **aminstarLmton;	// var to check message
   void aminstardecodeinit(double *y);
   int aminstardecode(unsigned char *c,unsigned long int maxnumloop,
							  unsigned long int &numloops);
   void aminstardecodeinit_old(double *y);
   int aminstardecode_old(unsigned char *c,unsigned long int maxnumloop,
							  unsigned long int &numloops);

   // variables and funtions for RCBP decoder
   void rcbpdecodeinit(double *y);
   int rcbpdecode(unsigned char *c,unsigned long int maxnumloop,
							  unsigned long int &numloops);
   
   // variables and funtions for bitflip1 decoder
   void bitflip1decodeinit(double *y,unsigned char *c);
   int bitflip1decode(unsigned char *c,unsigned long int maxnumloop,
							  unsigned long int &numloops);

   // variables and functions for GalAbitflip deceocer
   void galAbitflipdecodeinit(double *y);
   int galAbitflipdecode(unsigned char *c,
						 unsigned long int maxnumloop,
						 unsigned long int &numloops);
   unsigned char *yc;
   unsigned char **Bntom;
   unsigned char **Bmton;
   int *bmtonsum;
   int *bntomsum;

   // variables and functions for weighted bit flipping
   char *x;
   double *deltawbf;
   char *parities;
   double *betas;
   void weightedbitflipdecodeinit(double *y);
   int weightedbitflipdecode(double *y, unsigned char *c,
						 unsigned long int maxnumloop,
						 unsigned long int &numloops);
   

   // variable and functions for modified weighted bit flipping
   double modweightedbitflipalpha;
   void modweightedbitflipdecodeinit(double *y);
   int modweightedbitflipdecode(double *y, unsigned char *c,
						 unsigned long int maxnumloop,
						 unsigned long int &numloops);


   // gradient descent bit flip
   void grad1bitflipdecodeinit(double *y);
   int grad1bitflipdecode(double *y, unsigned char *c,
						 unsigned long int maxnumloop,
						 unsigned long int &numloops);

   // multibit gradient descent bit flip
   char *tryx;
   char *savex;
   void grad2bitflipdecodeinit(double *y);
   int grad2bitflipdecode(double *y, unsigned char *c,
						  unsigned long int maxnumloop,
						  unsigned long int &numloops);
   double gdobjectivefunct(double *y, char *tryx,
						   char *parities,
						   int &sumparities, double *deltawbf,
						   double &deltawbfmin, 
						   double &deltawbfmax, int &minidx);
   void gdcomputeparities(char *x, char *parities,
						  int &sumparities, bool &allparitiesgood);

   // Variables and functions for Divide and Concur
   double **DCrdat1;
   double **DCrdat2;
   double Lcsumsq;
   double EmaxDC;
   bool DCdoprobconstraint;
   
   void DC1decodeinit(double *y);
   int DC1decode(double *y, unsigned char *c, unsigned long int maxnumloop,
				unsigned long int &numloops);
   int DC1decodePD(double **rin, double **rout);
   void DC1decodePC(double **rin, double **rout, unsigned char *c);
   void DC1DM1(double **DCrdat1, double **DCrdat2);
   void DC1DM2(double **DCrdat1, double **DCrdat2);


   int grad2bitflip_mu;
   double grad2bitflip_theta;
   int *listtoflip;
   double fmin, fmax;
   double deltawbfmax, deltawbfmin;


   // DMBP
   double *DMBPb;				// beliefs on bits
   double DMBPZ;				// belief scale factor
   void DMBPdecodeinit(double *y);
   int DMBPdecode(unsigned char *c,
				  unsigned long int maxnumloop,
				  unsigned long int &numloops);

#ifdef DOLPSTUFF

   glp_prob *lp;    // linear programming structure
   glp_smcp glp_param; // linear programming parameters

   // Linear programming
   void gencomb(int n, int t, unsigned int *c, int *Nm,
				vector <vector <unsigned int> >&);
   // generate all combinations of t elements out of n
   unsigned int n_even_comb(unsigned int n,
							unsigned int& weightedncomb);

   unsigned int nchoosek(unsigned int n, unsigned int k);
   unsigned int factorial(unsigned int n);

   void LPinit(void);
   void LPdecodeinit(double *y);
   int LPdecode(unsigned char *c, unsigned long int maxnumloop,
				unsigned long int &numloops);


#endif // DOLPSTUFF
   

      // variables and functions for saving off some information
   int savestuff1;
   int savestuff2;
   std::ofstream fsave1;
   std::ofstream fsave2;

   void dosave1(const std::string fname1) {
	  savestuff1 = 1;
	  fsave1.open(fname1);
   };
   void dosave2(const std::string fname2) {
	  savestuff2 = 1;
	  fsave2.open(fname2);
   };
   void savestuff(double **Lmton, std::ofstream& fsave) {
	  int *na = new int[N];
	  bzero(na,N*sizeof(int));
	  for(int n = 0; n < N; n++) {
		 for(int m = 0; m < Mnlen[n]; m++) {
			fsave << Lmton[n][m] << std::endl;
		 }
	  }
	  delete[] na;
   }

};

#endif  // LDPCDECODER_H
/*
Local Variables:
compile-command: "g++ -c -g -std=c++11 ldpcdecoder.cc"
End:
*/

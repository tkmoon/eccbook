// ldpcdecoder.cc -- LDPC decoder class functions
// Todd K. Moon

#include "ldpcdecoder.h"
#include <fstream>
#include <iostream>
using namespace std;
#include <iomanip>
#include <chrono>
#include <stdlib.h>
#include <math.h>
#include <cfloat>
#include <algorithm> // (for min)
extern "C" {
#include <strings.h> // for bzero
}

int debugprint = 1;

// Adding a new decoder:
//  - create a #define for the name of the decoder in ldpcdecoder.h
//  - Add the name of the decoder to the variable decodernameslist[]
//         in ldpcdecoder.cc (below)
//  - Declare any class-level variables in ldpcdecoder.h
//  - Set the class-level array variables to 0 in LDPCDECODER::LDPCDECODER
//  - allocate space needed in LDPCDECODER::allocdecodedat() 
//    for the decoder
//  - free the space in LDPCDECODER::freedecodedat()
//  - Create an LDPCDECODER:: ...init() function (typically takes the
//     received data and sets up whatever is needed to
//     initialize the algorithm).  (e.g., LDPCDECODER::probdecodeinit( )
//     in ldpcdecoder.cc) and declare it in ldpcdecoder.h
//  - create a ...decode() function, e.g. LDPCDECODER::probdecode()
//     in ldpcdecoder.cc) and declare it in ldpcdecoder.h
// - in LDPCDECODER::decode(), add code to set the output array,
//     do timing, initialize, and decode, following the pattern there

std::string decodernameslist[] = {"Prob dec", // 0
								  "LL dec",   // 1
								  "Min sum",  // 2
								  "Phi decode",  // 3
								  "Q-Phi decode",  // 4
                                  "Min-Sum Correct",  // 5
								  "A-Min*",  // 6
								  "RCBP", // 7
								  "BitFlip1", // 8
								  "GalABitFlip", // 9
								  "WeightedBitFlip", // 10
								  "ModWeightedBitFlip", // 11
								  "GDBitFlip", // 12
								  "MultiGDbitflip",  // 13
								  "DivideConcur",  // 14
								  "Diff Map BP",  // 15
								  "Linear Programming", // 16
                                  };

int sigflag = 0; // used by signal catching function to break loops

LDPCDECODER::LDPCDECODER(string fname, int offset, unsigned long int indecodetype)
// fname = file containing sparse description of code
// offset = index offset in file.
//           0: C-like indexing   1: matlab-like indexing
// decodetype = DOPROBDEC for probability decoder,
//              DOLOGLIKEDECODE for log likelihood decoder
//              DOMINSUMDECODE for min sum decoder
//              DOPHIDECODE for phi decoder
//              DOQPHIDECODE for quantized phi decoder
//              DOMINSUMCORRECTDECODE min sum corrected decoder
//              DOAMINSTARDECODE approximate min* decoder
//              DOBCBPDECODE reduced-complexity box-plus
//              DOBITFLIP1 basic hard-decision bit flipping
//              GALABITFLIP basic hard-decision bit flipping
//              WEIGHTEDBITFLIP weighted bit flip
//              MODWEIGHTEDBITFLIP weighted bit flip
//              GDBITFLUP gradient descent bit flip
//              MULTIGDBITFLIP multi bit gradient descent bit flip
//              DC divide and concur   
//              DMBP difference-map belief propagation
//              LPDECODE linear programming decoding
// (OR multiple decodetypes together)
{
   int n,m,k,d;
   const int maxline = 1024;	// maximum length of line expected 
   char line[maxline];

   ifstream infile(fname);
   if(!infile) {
	  cerr << "Error: unable to open input file " << fname << endl;
	  exit(-1);
   }
   decodetype = indecodetype;

   // count the number of decoder types
   unsigned long int mask = 1;
   numdecodetype = 0;
   for(  ; mask <= decodetype;  ) {
	  if(mask & decodetype) numdecodetype++;
	  mask <<= 1;
   }

   // decodetypenum maps from the index 0.. numdecodetyp-1
   // to the decoder number.
   // decodernameslist[decodetypenum[i]] gives the name of the ith decoder

   // decodetypelist contains the bit map specification of the ith decoder

   decodetypelist = new int[numdecodetype];
   decodetypenum = new int[numdecodetype];
   for(int i=0, j = 0, mask = 1; mask <= decodetype; mask <<= 1) {
	  if(mask & decodetype) {
		 decodetypelist[i] = mask;
		 decodetypenum[i] = j;
		 i++;
	  }
	  j++;
   }

   infile >> N;
   infile >> M;
   K = N-M;
   infile >> maxcolwt;
   infile >> maxrowwt;

   CALLOCMATRIX(Mn,int,N,maxcolwt);
   CALLOCMATRIX(Nm,int,M,maxrowwt);

   // read in the column weights
   Mnlen = new int[N];
   for(n = 0; n < N; n++) {
	  infile >> Mnlen[n];
   }
   // read in the row weights
   Nmlen = new int[M];
   for(m = 0; m < M; m++) {
	  infile >> Nmlen[m];
   }
   // read in the Mn data (checks for each bit)
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 infile >> d;
		 Mn[n][m] = d - offset;
	  }
   }

   // read in the Nm data (bits for each check)
   for(m = 0; m < M; m++) {
	  for(n = 0; n < Nmlen[m]; n++) {
		 infile >> d;
		 Nm[m][n] = d-offset;
	  }
   }
   allocdecodedat();
   if(decodetype & DOQPHIDECODE) {
	  buildqphitable();		// build the quantized phi table
   }
   savestuff1 = 0;  savestuff2 = 0;
}


void
LDPCDECODER::allocdecodedat(void)
// allocate memory used in the decoder
{
   int i,j;

   CALLOCMATRIX(c, unsigned char,numdecodetype,N); // decoded output
   problikedata = new double[N];  // probability or likelihood data
   na = new unsigned int[N];  // "number above" data
   int nalloc = MAX(maxcolwt, maxrowwt); 
   // data pulled from rows to use with forward/backward tables
   rowdata = new double[maxrowwt];
   // data used for forward/backward tables

   forwardprod = new double[nalloc-1];  // used for/back algs for all decoders
   backwardprod = new double[nalloc-1];
   forbackdataout = new double[nalloc];

   outproblikedata1 = new double[N];  // output probability or likelhood data
   outproblikedata0 = 0; // used by prob -- may not be needed

   CALLOCMATRIX(mtondata1,double,N,maxcolwt);  // m to n message data
   CALLOCMATRIX(ntomdata1,double,maxcolwt,N);  // n to m message data
   mtondata0 = 0;  // used by prob -- may not be needed 
   ntomdata0 = 0;  // usd by prob -- may not be needed

   Pntom0 = 0;					// used by prob
   deltar = 0;					// used by prob
   Pnout0 = 0;					// used by prob
   poutprods1 = 0;				// used by prob
   forwardprod1 = 0;			// used by prob
   backwardprod1 = 0;			// used by prob

   alphantom = 0;				// used by phi-based decoders
   alpharow = 0;				// used by phi-based decoders

   synd = 0;					   // used by bit flipping
   flist = 0;					   // list of syndromes with largest count


   yc = 0;						// used by Gallager Algorithm A
   Bntom = 0;
   Bmton = 0;
   bntomsum = 0;

   x = 0;						// used by weighted bit flip
   parities = 0;
   betas = 0;
   deltawbf = 0;
   tryx = 0;
   savex = 0;
   
   listtoflip = 0;

   DCrdat1 = 0;
   DCrdat2 = 0;
   Lc = 0;
   // DCdoprobconstraint = false;
   DCdoprobconstraint = true;

   if(decodetype & DOPROBDECODE) {			// probability decoder
	  pn = problikedata;
	  deltaq = rowdata; // new double[maxrowwt];

	  Pmton1 = mtondata1;   // N, maxcolwt
	  CALLOCMATRIX(Pmton0,double,N,maxcolwt);

	  Pntom1 = ntomdata1;  // pmaxcolwt, N
	  CALLOCMATRIX(Pntom0,double,maxcolwt,N);

	  CALLOCMATRIX(deltar,double,maxcolwt,N); // used only for prob decoder
	  Pnout1 = outproblikedata1;
	  Pnout0 = new double[N];
	  
	  poutprods = forbackdataout;   // used for CN to VN step

	  poutprods0 = forbackdataout;  // (re) used for VN to CN step
	  poutprods1 = new double[nalloc];
	  
	  forwardprod0 = forwardprod;  // reference memory already allocated
	  backwardprod0 = backwardprod;

	  forwardprod1 = new double[maxrowwt];
	  backwardprod1 = new double[maxrowwt];
   }
   if(decodetype & DOLOGLIKEDECODE) { // loglikelhood decoder
	  Lc = problikedata;
	  tanhstuff = rowdata;
	  poutprobsLL = forbackdataout;
	  Lcout = outproblikedata1;

	  Lmton = mtondata1;
	  Lntom = ntomdata1;
	  poutsum = forbackdataout;
	  forwardsum0 = forwardprod;   // reference memory already allocated
	  backwardsum0 = backwardprod;
   }
   if(decodetype & DOMINSUMDECODE) { // minsum decoder
	  if(!alphantom) {CALLOCMATRIX(alphantom,char,maxcolwt,N); }
	  if(!alpharow) alpharow = new char[maxrowwt];

	  MSLc = problikedata;
	  MSLcout = outproblikedata1;
	  betantom = ntomdata1;
	  MSLmton = mtondata1;
	  
	  betarow = rowdata;
	  poutmin = forbackdataout;
	  poutsum = forbackdataout;
	  forwardmin0 = forwardprod;   // reference memory already allocated
	  backwardmin0 = backwardprod;
	  forwardsum0 = forwardprod;   // reference memory already allocated
	  backwardsum0 = backwardprod;
	  poutsum = forbackdataout;  // new double[maxrowwt];
   }
   if(decodetype & (DOPHIDECODE | DOQPHIDECODE)) { // phi or quantized phi
	  if(!alphantom) { CALLOCMATRIX(alphantom,char,maxcolwt,N);}
	  if(!alpharow) alpharow = new char[maxrowwt];
	  phiLc = problikedata;
	  phiLcout = outproblikedata1;
	  betantom = ntomdata1;
	  phiLmton = mtondata1;
	  
	  betarow = rowdata;
	  poutmin = forbackdataout;
	  poutsum = forbackdataout;
	  forwardmin0 = forwardprod;   // reference memory already allocated
	  backwardmin0 = backwardprod;
	  forwardsum0 = forwardprod;   // reference memory already allocated
	  backwardsum0 = backwardprod;
	  poutsum = forbackdataout;  // new double[maxrowwt];
	  phibetarow = rowdata;
	  phioutphisum = forbackdataout;
   }
   if(decodetype & DOMINSUMCORRECTDECODE) { // minsum decoder
	  if(!alphantom) {CALLOCMATRIX(alphantom,char,maxcolwt,N); }
	  if(!alpharow) alpharow = new char[maxrowwt];

	  MScorrectLc = problikedata;
	  MScorrectLcout = outproblikedata1;
	  betantom = ntomdata1;
	  MScorrectLmton = mtondata1;
	  
	  betarow = rowdata;
	  poutmin = forbackdataout;
	  poutsum = forbackdataout;
	  forwardmin0 = forwardprod;   // reference memory already allocated
	  backwardmin0 = backwardprod;
	  forwardsum0 = forwardprod;   // reference memory already allocated
	  backwardsum0 = backwardprod;
	  poutsum = forbackdataout;  // new double[maxrowwt];
   }
   if(decodetype & DOAMINSTARDECODE) { // minsum decoder
	  if(!alphantom) {CALLOCMATRIX(alphantom,char,maxcolwt,N); }
	  if(!alpharow) alpharow = new char[maxrowwt];

	  aminstarLc = problikedata;
	  aminstarLcout = outproblikedata1;
	  betantom = ntomdata1;
	  aminstarLmton = mtondata1;
	  
	  betarow = rowdata;
	  poutmin = forbackdataout;
	  poutsum = forbackdataout;
	  forwardmin0 = forwardprod;   // reference memory already allocated
	  backwardmin0 = backwardprod;
	  forwardsum0 = forwardprod;   // reference memory already allocated
	  backwardsum0 = backwardprod;
	  poutsum = forbackdataout;  // new double[maxrowwt];
   }
   if(decodetype & DORCBPDECODE) { // RCBP decoder
	  if(!alphantom) {CALLOCMATRIX(alphantom,char,maxcolwt,N); }
	  if(!alpharow) alpharow = new char[maxrowwt];

	  aminstarLc = problikedata;
	  aminstarLcout = outproblikedata1;
	  betantom = ntomdata1;
	  aminstarLmton = mtondata1;
	  
	  betarow = rowdata;
	  poutmin = forbackdataout;
	  poutsum = forbackdataout;
	  forwardmin0 = forwardprod;   // reference memory already allocated
	  backwardmin0 = backwardprod;
	  forwardsum0 = forwardprod;   // reference memory already allocated
	  backwardsum0 = backwardprod;
	  poutsum = forbackdataout;  // new double[maxrowwt];
   }
   if(decodetype & DOBITFLIP1) {
	  if(!synd) { synd = new unsigned char[M]; }
	  if(!flist) {flist = new unsigned int[N]; }  // way more than needed
   }
   if(decodetype & GALABITFLIP) {
	  if(!yc) yc = new unsigned char[N];
	  if(!Bntom) { CALLOCMATRIX(Bntom,unsigned char,maxcolwt,N);}
	  if(!Bmton) { CALLOCMATRIX(Bmton,unsigned char,maxcolwt,N);}
	  if(!bntomsum) bntomsum = new int[N];
   }
   if(decodetype & (WEIGHTEDBITFLIP | MODWEIGHTEDBITFLIP)) {
	  if(!x) x = new char[N];
	  if(!deltawbf) deltawbf = new double[N];
	  if(!parities) parities = new char[M];
	  if(!betas) betas = new double[M];
   }
   if(decodetype & GDBITFLIP) {
	  if(!deltawbf) deltawbf = new double[N];
	  if(!x) x = new char[N];
	  if(!parities) parities = new char[M];
   }
   if(decodetype & MULTIGDBITFLIP) {
	  if(!deltawbf) deltawbf = new double[N];
	  if(!x) x = new char[N];
	  if(!tryx) tryx = new char[N];
	  if(!savex) savex = new char[N];
	  if(!parities) parities = new char[M];
	  if(!listtoflip) listtoflip = new int[N];
   }

   if(decodetype & DC1) {  // divide and concur
	  if(!DCrdat1) { CALLOCMATRIX(DCrdat1,double,(maxcolwt+1),N);}
	  if(!DCrdat2) { CALLOCMATRIX(DCrdat2,double,(maxcolwt+1),N);}
	  if(!Lc)  Lc  = new double[N];
	  if(!x) x = new char[N];
   }

   if(decodetype & DMBP) {
	  if(!alphantom) {CALLOCMATRIX(alphantom,char,maxcolwt,N); }
	  if(!alpharow) alpharow = new char[maxrowwt];
	  betarow = rowdata;
	  betantom = ntomdata1;
	  MSLc = problikedata;
	  DMBPb = outproblikedata1;
	  poutmin = forbackdataout;
	  betantom = ntomdata1;
	  MScorrectLmton = mtondata1;
	  forwardmin0 = forwardprod;   // reference memory already allocated
	  backwardmin0 = backwardprod;
	  DMBPZ = 0.445;
   }


#ifdef DOLPSTUFF
   if(decodetype & LPDECODE) {
	  LPinit();
   }

#endif
}

void
LDPCDECODER::freedecodedat(void)
{
   FREEMATRIX(c);
   delete [] problikedata;
   delete [] na;
   delete [] rowdata;
   delete []forwardprod;
   delete [] backwardprod;
   delete []outproblikedata1;
   FREEMATRIX(mtondata1);
   FREEMATRIX(ntomdata1);
   if(Pntom0) {FREEMATRIX(Pntom0); }
   if(deltar) {FREEMATRIX(deltar); }
   if(Pnout0) delete[] Pnout0;
   if(poutprods1) delete[] poutprods1;
   if(forwardprod1) delete[] forwardprod1;
   if(backwardprod1) delete[] backwardprod1;
   if(alphantom) {FREEMATRIX(alphantom)};
   if(alpharow) delete[] alpharow;
   if(synd) delete[] synd;
   if(flist) delete[] flist;
   if(yc) delete[] yc;
   if(Bntom) { FREEMATRIX(Bntom); }
   if(Bmton) { FREEMATRIX(Bmton); }
   if(bntomsum) delete[] bntomsum;
   
   if(x) delete[] x;
   if(deltawbf) delete[] deltawbf;
   if(parities) delete[] parities;
   if(betas) delete[] betas;
   if(listtoflip) delete[] listtoflip;
   if(tryx) delete[] tryx;
   if(savex) delete[] savex;
   
   if(DCrdat1) {FREEMATRIX(DCrdat1)};
   if(DCrdat2) {FREEMATRIX(DCrdat2)};
}

int LDPCDECODER::decodetypetonum(unsigned long int decodetype)
{
   for(int i = 0; i < numdecodetype; i++) {
	  if(decodetypelist[i] == decodetype) {
		 return i;
	  }
   }
   return -1;

}
   
unsigned long int
LDPCDECODER::decode(double *y, int inprintstuff, unsigned long int maxnumloop,
					unsigned long int *numloops,double *durations)
// y = output of channel 
// printstuff - set to print intermediate results
// maxnumloop = maximum number of decoding iterations
// numloops (one for each decoder): indicat number of loops
//   if numloops[.] is a positive integer on input, then
//   the decoder does that number of iterations.
//   if numloops[.] is 0 in input, then the deocder proceeds
//   until parities all check or until maxnumloop iterations
// durations is duration (in seconds) of decoding (for each decoder type)

// returns: 1 in bit position for decoder type if decoding succeeds;
// 0 for a decoding failure
{
   unsigned long int paritychecks = 0;
   unsigned char *cc;			// pointer to decoded bits for this decoder
   int parity;
   int decodenum;
   
   sigflag = 0;
   printstuff = inprintstuff;

   if(decodetype & DOPROBDECODE) {
	  decodenum = decodetypetonum(DOPROBDECODE);
	  probdecodeinit(y);
	  cc = c[decodenum];
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = probdecode(cc,maxnumloop,
							   numloops[decodenum]);
//cout << "!" << flush;
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];
	  if(parity) { // correctly decoded
		 paritychecks |= DOPROBDECODE;
	  }
   }

   if(decodetype & DOLOGLIKEDECODE) {
	  decodenum = decodetypetonum(DOLOGLIKEDECODE);
	  lldecodeinit(y);
	  cc = c[decodenum];
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = lldecode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];
	  if(parity) { // correctly decoded
		 paritychecks |= DOLOGLIKEDECODE;
	  }
   }

   if(decodetype & DOMINSUMDECODE) {
	  decodenum = decodetypetonum(DOMINSUMDECODE);
	  minsumdecodeinit(y);
	  cc = c[decodenum];
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = minsumdecode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= DOMINSUMDECODE;
	  }
   }

   if(decodetype & DOPHIDECODE) {
	  decodenum = decodetypetonum(DOPHIDECODE);
	  phidecodeinit(y);
	  cc = c[decodenum];
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = phidecode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= DOPHIDECODE;
	  }
   }

   if(decodetype & DOQPHIDECODE) {
	  decodenum = decodetypetonum(DOQPHIDECODE);
	  phidecodeinit(y);   // same for phi and qphi decoders
	  cc = c[decodenum];
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = qphidecode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= DOQPHIDECODE;
	  }
   }

   if(decodetype & DOMINSUMCORRECTDECODE) {
	  decodenum = decodetypetonum(DOMINSUMCORRECTDECODE);
	  minsumdecodeinit(y);
	  cc = c[decodenum];
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = minsumcorrectdecode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= DOMINSUMCORRECTDECODE;
	  }
   }

   if(decodetype & DOAMINSTARDECODE) {
	  decodenum = decodetypetonum(DOAMINSTARDECODE);
	  aminstardecodeinit(y);
	  cc = c[decodenum];
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = aminstardecode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= DOAMINSTARDECODE;
	  }
   }

   if(decodetype & DORCBPDECODE) {
	  decodenum = decodetypetonum(DORCBPDECODE);
	  rcbpdecodeinit(y);
	  cc = c[decodenum];
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = rcbpdecode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= DORCBPDECODE;
	  }
   }

   if(decodetype & DOBITFLIP1) {
	  decodenum = decodetypetonum(DOBITFLIP1);
	  cc = c[decodenum];
	  bitflip1decodeinit(y,cc);
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = bitflip1decode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= DOBITFLIP1;
	  }
   }
   if(decodetype & GALABITFLIP) {
	  decodenum = decodetypetonum(GALABITFLIP);
	  cc = c[decodenum];
	  galAbitflipdecodeinit(y);
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = galAbitflipdecode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= GALABITFLIP;
	  }
   }

   if(decodetype & WEIGHTEDBITFLIP) {
	  decodenum = decodetypetonum(WEIGHTEDBITFLIP);
	  cc = c[decodenum];
	  weightedbitflipdecodeinit(y);
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = weightedbitflipdecode(y,cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= WEIGHTEDBITFLIP;
	  }
   }

   if(decodetype & MODWEIGHTEDBITFLIP) {
	  decodenum = decodetypetonum(MODWEIGHTEDBITFLIP);
	  cc = c[decodenum];
	  modweightedbitflipdecodeinit(y);
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = modwtedbitflipdecode(y,cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= MODWEIGHTEDBITFLIP;
	  }
   }

   if(decodetype & GDBITFLIP) {
	  decodenum = decodetypetonum(GDBITFLIP);
	  cc = c[decodenum];
	  grad1bitflipdecodeinit(y);
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = grad1bitflipdecode(y,cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= GDBITFLIP;
	  }
   }

   if(decodetype & MULTIGDBITFLIP) {
	  decodenum = decodetypetonum(MULTIGDBITFLIP);
	  cc = c[decodenum];
	  grad2bitflipdecodeinit(y);
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = grad2bitflipdecode(y,cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= MULTIGDBITFLIP;
	  }
   }

   if(decodetype & DC1) {
	  decodenum = decodetypetonum(DC1);
	  cc = c[decodenum];
	  DC1decodeinit(y);
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = DC1decode(y,cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= DC1;
	  }
   }

   if(decodetype & DMBP) {
	  decodenum = decodetypetonum(DMBP);
	  cc = c[decodenum];
	  DMBPdecodeinit(y);
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = DMBPdecode(cc,maxnumloop,numloops[decodenum]);
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration)/numloops[decodenum];

	  if(parity) { // correctly decoded
		 paritychecks |= DMBP;
	  }
   }

#ifdef DOLPSTUFF
   if(decodetype & LPDECODE) {
	  decodenum = decodetypetonum(LPDECODE);
	  cc = c[decodenum];

	  LPdecodeinit(y);
	  chrono::high_resolution_clock::time_point t1 =
		 chrono::high_resolution_clock::now();
	  parity = LPdecode(cc,maxnumloop,numloops[decodenum]);
//cout << "." << flush;
	  chrono::high_resolution_clock::time_point t2 =
		 chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>
		 ( t2 - t1 ).count();
	  durations[decodenum] += double(duration);  // /numloops[decodenum];
	  // Note:  Since there are not iterations in the usual sense,
	  // duration is not normalized by number of loops.

	  if(parity) { // correctly decoded
		 paritychecks |= LPDECODE;
	  }
   }

#endif
   
   return(paritychecks);
}


int LDPCDECODER::checkparity(unsigned char *c)
// c = (alleged) codeword
// returns 1 if all the parities check
// returns 0 if any of the parities do not check
{
   int parity = 1;
   int m,l;
   unsigned int z;
   
   parity = 1;
   for(m = 0; m < M; m++) { // check all parities
	  z = 0;
	  for(l = 0; l < Nmlen[m]; l++) {
		 z += c[Nm[m][l]];
	  }
	  z %= 2;
	  if(z) {
		 // Parity check fails.  Bail out of check at this point.
		 parity = 0;
		 break;
	  }
   }
   return parity;
}   

void LDPCDECODER::probdecodeinit(double *y)
   // convert to posterior probabilities
{
   int i,m,n;
   if(printstuff) VECDUMP(y,N);
   double L;

   L = 2*a/sigma2;
   // compute channel posterior probabilities
   for(i = 0; i < N; i++) {
	  pn[i] = 1/(1+exp(L*y[i]));
   }

   if(printstuff) VECDUMP(pn,N);
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 Pntom1[m][n] = pn[n];
	  }
   }
   if(printstuff) printsparse("q(init)",Pntom1);
}



int LDPCDECODER::probdecode(unsigned char *c, unsigned long int maxnumloop,
							unsigned long int &numloops)
// returns 1 if parity checks, 0 if parity does not check
{
   int l,m,n;
   double pall0, pall1;
   double tp;
   int parity = 0;
   unsigned long int loopcount = 0;

   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do { // loop until completion of decoding, or loop count
	  loopcount++;
	  // initialize the indexing for dealing with sparse matrices

	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // Check node to Variable node step
	  for(m = 0; m < M; m++) {  // for each check node
		 // compute the delta messages
		 for(l = 0; l < Nmlen[m]; l++) {
			n = Nm[m][l];	  // column index of lth element on this row
			deltaq[l] = 1-2*Pntom1[na[n]][n]; // compute delta q
		 }

		 // compute leave-one-out products on delta messages
		 forbackprobs(deltaq,Nmlen[m], poutprods);
		 for(l = 0; l < Nmlen[m]; l++) {  // for each variable node on this row
			n = Nm[m][l];

			// normalize and assign back into sparse structure
			Pmton1[n][na[n]] = (1-poutprods[l])/2;   //   (1-prod)/2;
			Pmton0[n][na[n]++] = (1+poutprods[l])/2;
		 }
	  } // for m
	  if(printstuff) printsparsetranspose("Pmton1",Pmton1);
 
	  // Variable node to Check node step
	  for(n = 0; n < N; n++) {  // for each variable node
		 // compute leave-one-out product on Pmton0 and Pmton1, including
		 // 1-pn[n] and pn[n] as initial factors
		 forbackprobs2(1-pn[n],Pmton0[n],pn[n],Pmton1[n],Mnlen[n],poutprods0,
					   pall0,poutprods1,pall1);
		 // normalize and assign output probabilities
		 tp = pall0 + pall1;
		 Pnout0[n] = pall0/tp;   // produce normalized output probabilities
		 Pnout1[n] = pall1/tp;
		 // normalize and assign each message to variable node
		 for(l = 0; l < Mnlen[n]; l++) {
			tp = poutprods0[l] + poutprods1[l];
			Pntom0[l][n] = poutprods0[l]/tp;
			Pntom1[l][n] = poutprods1[l]/tp;
		 }
		 // threshold output probabilities to get decoded bits
		 if(Pnout1[n] > 0.5) c[n] = 1; else c[n] = 0;
	  } // for n
	  
	  if(printstuff) printsparse("Pntom1",Pntom1);
	  if(printstuff) VECDUMP(Pnout1,N);
	  if(printstuff) VECDUMP2(c,N,int);

	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
			if(parity) { // if parities check
			   numloops = loopcount;
			   return parity;
			}
	  }
	  if(sigflag) {
		 cout << "sigflag set: breaking" << endl;
		 break;
	  }
   } while(loopcount < maxnumloop);
   if(numloops) {   // parity has not been previously checked
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}


void LDPCDECODER::lldecodeinit(double *y)
{
   int n,m;
   double Lc1;
   Lc1 = 2*a/sigma2;// a is negative if the bit 0 is represented with a negative

   for(n = 0; n < N; n++) {
	  Lc[n] = Lc1*y[n];
   }
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 Lntom[m][n] = Lc[n];
	  }
   }
   if(printstuff) printsparse("Lctom(init)",Lntom);
}

int LDPCDECODER::lldecode(unsigned char *c, unsigned long int maxnumloop,
							unsigned long int &numloops)

{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;

   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check
		 // work over nonzero elements of this row
		 for(l = 0; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			tanhstuff[l] = tanh(Lntom[na[n]][n]/2);
		 }
		 // compute leave-one-out products on tanh stuff
		 forbackprobs(tanhstuff,Nmlen[m],poutprobsLL);
		 for(l = 0; l < Nmlen[m]; l++) { // for each n on this row
			n = Nm[m][l];
			Lmton[n][na[n]++] = 2*atanh(poutprobsLL[l]);
		 }
	  }  // for m
	  
      if(printstuff) printsparsetranspose("Lmton",Lmton);
	  if(savestuff1) savestuff(Lmton,fsave1);
   
	  // Variable node to check node computation
	  for(n = 0; n < N; n++) {  // work over each column
		 forbacksum(Lc[n],Lmton[n],Mnlen[n],poutsum,sumall);
		 Lcout[n] = sumall;
		 if(Lcout[n] >= 0) c[n] = 0; else c[n] = 1;

		 for(l = 0; l < Mnlen[n]; l++) {
			Lntom[l][n] = poutsum[l];
		 }
	  } // for n
	  if(printstuff) printsparse("Lntom",Lntom);
	  if(printstuff) VECDUMP(Lcout,N);
	  if(printstuff) VECDUMP2(c,N,int);

	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
		 if(parity) { // if parities check
			numloops = loopcount;
			return parity;
		 }
	  }
   } while(loopcount < maxnumloop);
   if(numloops) { // parity has not been previously checked
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}


void LDPCDECODER::minsumdecodeinit(double *y)
{
   int n,m;
   double Lc;
   Lc = 2*a/sigma2;
   for(n = 0; n < N; n++) {
	  MScorrectLc[n] = Lc*y[n];
   }
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 alphantom[m][n] = (MScorrectLc[n] < 0);
		 betantom[m][n] = abs(MScorrectLc[n]);
	  }
   }
   if(printstuff) printsparse("MS Lntom(beta,init)",betantom);
   if(printstuff) printsparsechar("alphantom",alphantom);

}


int LDPCDECODER::minsumdecode(unsigned char *c, unsigned long int maxnumloop,
							unsigned long int &numloops)

{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;
   char alphasign, alphacum;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check
		 // work over nonzero elements of this row
		 l = 0;
		 n = Nm[m][l];
		 alphasign = alphantom[na[n]][n];
		 alpharow[l] = alphasign;
		 alphacum = alphasign;
		 betarow[l] = betantom[na[n]][n];
		 for(l = 1; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			alphasign = alphantom[na[n]][n];
			alphacum ^= alphasign;
			alpharow[l] = alphasign;
			betarow[l] = betantom[na[n]][n];
		 }
		 // compute leave-one-out min data
		 forbackmin(betarow,Nmlen[m],poutmin);
		 for(l = 0; l < Nmlen[m]; l++) { // for each n on this row
			n = Nm[m][l];
			alphasign = alphacum ^ alphantom[na[n]][n]; 
			MSLmton[n][na[n]++] = (1-2*alphasign)*poutmin[l];
		 }
	  }  // for m
	  
      if(printstuff) printsparsetranspose("MSLmton",MSLmton);
	  if(savestuff1) savestuff(MSLmton,fsave2);

	  // Variable node to check node computation
	  for(n = 0; n < N; n++) {  // work over each column
		 forbacksum(MSLc[n],MSLmton[n],Mnlen[n],poutsum,sumall);
		 MSLcout[n] = sumall;
		 if(MSLcout[n] >= 0) c[n] = 0; else c[n] = 1;
		 for(l = 0; l < Mnlen[n]; l++) {
			betantom[l][n] = abs(poutsum[l]);
			alphantom[l][n] = poutsum[l] < 0;
		 }
	  } // for n
	  if(printstuff) printsparse("betantom",betantom);
	  if(printstuff) printsparsechar("alphantom",alphantom);
	  if(printstuff) VECDUMP(MSLcout,N);
	  if(printstuff) VECDUMP2(c,N,int);

	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
		 if(parity) { // if parities check
			numloops = loopcount;
			return parity;
		 }
	  }
   } while(loopcount < maxnumloop);
   if(numloops) { // parity has not been previously checked
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}


void LDPCDECODER::minsumcorrectdecodeinit(double *y)
{
   int n,m;
   double Lc;
   Lc = 2*a/sigma2;
   for(n = 0; n < N; n++) {
	  MSLc[n] = Lc*y[n];
   }
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 alphantom[m][n] = (MScorrectLc[n] < 0);
		 betantom[m][n] = abs(MScorrectLc[n]);
	  }
   }
   if(printstuff) printsparse("MS correct Lntom(beta,init)",betantom);
   if(printstuff) printsparsechar("alphantom",alphantom);

}


int LDPCDECODER::minsumcorrectdecode(unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)

{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;
   char alphasign, alphacum;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check
		 // work over nonzero elements of this row
		 l = 0;
		 n = Nm[m][l];
		 alphasign = alphantom[na[n]][n];
		 alpharow[l] = alphasign;
		 alphacum = alphasign;
		 betarow[l] = betantom[na[n]][n];
		 for(l = 1; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			alphasign = alphantom[na[n]][n];
			alphacum ^= alphasign;
			alpharow[l] = alphasign;
			betarow[l] = betantom[na[n]][n];
		 }
		 // compute leave-one-out min data
		 forbackmincorrect(betarow,Nmlen[m],poutmin);
		 for(l = 0; l < Nmlen[m]; l++) { // for each n on this row
			n = Nm[m][l];
			alphasign = alphacum ^ alphantom[na[n]][n]; 
			MScorrectLmton[n][na[n]++] = (1-2*alphasign)*poutmin[l];
		 }
	  }  // for m
	  
      if(printstuff) printsparsetranspose("MScorrectLmton",MScorrectLmton);
	  if(savestuff1) savestuff(MScorrectLmton,fsave2);

	  // Variable node to check node computation
	  for(n = 0; n < N; n++) {  // work over each column
		 forbacksum(MScorrectLc[n],MScorrectLmton[n],Mnlen[n],poutsum,sumall);
		 MScorrectLcout[n] = sumall;
		 if(MScorrectLcout[n] >= 0) c[n] = 0; else c[n] = 1;
		 for(l = 0; l < Mnlen[n]; l++) {
			betantom[l][n] = abs(poutsum[l]);
			alphantom[l][n] = poutsum[l] < 0;
		 }
	  } // for n
	  if(printstuff) printsparse("betantom",betantom);
	  if(printstuff) printsparsechar("alphantom",alphantom);
	  if(printstuff) VECDUMP(MScorrectLcout,N);
	  if(printstuff) VECDUMP2(c,N,int);

	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
		 if(parity) { // if parities check
			numloops = loopcount;
			return parity;
		 }
	  }
   } while(loopcount < maxnumloop);
   if(numloops) { // parity has not been previously checked
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}

void LDPCDECODER::aminstardecodeinit_old(double *y)
{
   int n,m;
   double Lc;
   Lc = 2*a/sigma2;
   for(n = 0; n < N; n++) {
	  aminstarLc[n] = Lc*y[n];
   }
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 alphantom[m][n] = (aminstarLc[n] < 0);
		 betantom[m][n] = abs(aminstarLc[n]);
	  }
   }
   if(printstuff) printsparse("aminstarLntom(beta,init)",betantom);
   if(printstuff) printsparsechar("alphantom",alphantom);

}


int LDPCDECODER::aminstardecode_old(unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)

{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;
   char alphasign, alphacum;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array
	  int nmin;							  // index with smallest reliability
	  
	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check
		 double beta1min = DBL_MAX;
		 int lmin;
		 double sumphi, savephi, phitotal, phisumphi;
		 // work over nonzero elements of this row
		 l = 0;
		 n = Nm[m][l];
		 alphasign = alphantom[na[n]][n];
		 alpharow[l] = alphasign;
		 alphacum = alphasign;
		 betarow[l] = betantom[na[n]][n];
		 for(l = 1; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			alphasign = alphantom[na[n]][n];
			alphacum ^= alphasign;
			alpharow[l] = alphasign;
			betarow[l] = betantom[na[n]][n];
			if(betarow[l] < beta1min) {
			   beta1min = betarow[l];
			   nmin = n;
			   lmin = l;
			}
		 }
		 sumphi = 0;
		 for(l = 0; l < Nmlen[m]; l++) {
			if(l == lmin) {
			   savephi = phi(betarow[l]);
			}
			else {
			   sumphi += phi(betarow[l]);
			}
		 }
		 phitotal = fabs(phi(savephi + sumphi));
		 phisumphi = fabs(phi(sumphi));
		 // the two messages are phi(sumphi) and (sumphi + savephi)
		 for(l = 0; l < Nmlen[m]; l++) {  // message to each V_n
			n = Nm[m][l];
			alphasign = alphacum ^ alphantom[na[n]][n];
			if(l == lmin) {
			   aminstarLmton[n][na[n]++] = (1 - 2*alphasign)*phisumphi;
			}
			else {
			   aminstarLmton[n][na[n]++] = (1 - 2*alphasign)*phitotal;
			}
		 }
	  }
      if(printstuff) printsparsetranspose("aminstarcorrectLmton",aminstarLmton);

	  // Variable node to check node computation
	  for(n = 0; n < N; n++) {  // work over each column
		 forbacksum(aminstarLc[n],aminstarLmton[n],Mnlen[n],poutsum,sumall);
		 aminstarLcout[n] = sumall;
		 if(aminstarLcout[n] >= 0) c[n] = 0; else c[n] = 1;
		 for(l = 0; l < Mnlen[n]; l++) {
			betantom[l][n] = abs(poutsum[l]);
			alphantom[l][n] = poutsum[l] < 0;
		 }
	  } // for n
	  if(printstuff) printsparse("betantom",betantom);
	  if(printstuff) printsparsechar("alphantom",alphantom);
	  if(printstuff) VECDUMP(MScorrectLcout,N);
	  if(printstuff) VECDUMP2(c,N,int);

	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
		 if(parity) { // if parities check
			numloops = loopcount;
			return parity;
		 }
	  }
   } while(loopcount < maxnumloop);
   if(numloops) { // parity has not been previously checked
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}


void LDPCDECODER::aminstardecodeinit(double *y)
{
   int n,m;
   double Lc;
   Lc = 2*a/sigma2;
   for(n = 0; n < N; n++) {
	  aminstarLc[n] = Lc*y[n];
   }
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 alphantom[m][n] = (aminstarLc[n] < 0);
		 betantom[m][n] = abs(aminstarLc[n]);
	  }
   }
   if(printstuff) printsparse("aminstarLntom(beta,init)",betantom);
   if(printstuff) printsparsechar("alphantom",alphantom);

}

int LDPCDECODER::aminstardecode(unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)

{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;
   char alphasign, alphacum;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array
	  int nmin;							  // index with smallest reliability
	  
	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check
		 double Lnmin2m = DBL_MAX;
		 double Lm2nmin = DBL_MAX;
		 double Lm2star;
		 int alphacum = 0;
		 for(l = 0; l < Nmlen[m]; l++) { // for each adjacent variable node
			n = Nm[m][l];
			alphacum = alphacum ^ alphantom[na[n]][n];
			if(betantom[na[n]][n] < Lnmin2m) { // new smallest incoming
			   Lm2nmin = boxplus(Lm2nmin,Lnmin2m);
			   nmin = n;
			   Lnmin2m = betantom[na[n]][n];
			}
			else { // not a new min
			   Lm2nmin = boxplus(Lm2nmin,betantom[na[n]][n]);
			}
		 } // for
		 Lm2star = boxplus(Lm2nmin,Lnmin2m);

		 // compute the outgoing messages
		 for(l = 0; l < Nmlen[m]; l++) { // for each adjacent variable node
			n = Nm[m][l];
			alphasign = alphacum ^ alphantom[na[n]][n];
			if(n == nmin) { // message to least reliable
			   aminstarLmton[n][na[n]++] = (1 - 2*alphasign)*Lm2nmin;
			}
			else {
			   aminstarLmton[n][na[n]++] = (1 - 2*alphasign)*Lm2star;
			}
		 }
	  } // for m
		 

      if(printstuff) printsparsetranspose("aminstarcorrectLmton",aminstarLmton);

	  // Variable node to check node computation
	  for(n = 0; n < N; n++) {  // work over each column
		 forbacksum(aminstarLc[n],aminstarLmton[n],Mnlen[n],poutsum,sumall);
		 aminstarLcout[n] = sumall;
		 if(aminstarLcout[n] >= 0) c[n] = 0; else c[n] = 1;
		 for(l = 0; l < Mnlen[n]; l++) {
			betantom[l][n] = abs(poutsum[l]);
			alphantom[l][n] = poutsum[l] < 0;
		 }
	  } // for n
	  if(printstuff) printsparse("betantom",betantom);
	  if(printstuff) printsparsechar("alphantom",alphantom);
	  if(printstuff) VECDUMP(MScorrectLcout,N);
	  if(printstuff) VECDUMP2(c,N,int);

	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
		 if(parity) { // if parities check
			numloops = loopcount;
			return parity;
		 }
	  }
   } while(loopcount < maxnumloop);
   if(numloops) { // parity has not been previously checked
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}



void LDPCDECODER::rcbpdecodeinit(double *y)
{
   int n,m;
   double Lc;
   Lc = 2*a/sigma2;
   for(n = 0; n < N; n++) {
	  aminstarLc[n] = Lc*y[n];
   }
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 alphantom[m][n] = (aminstarLc[n] < 0);
		 betantom[m][n] = abs(aminstarLc[n]);
	  }
   }
   if(printstuff) printsparse("aminstarLntom(beta,init)",betantom);
   if(printstuff) printsparsechar("alphantom",alphantom);
}



int LDPCDECODER::rcbpdecode(unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)
{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;
   char alphasign, alphacum;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array
	  int nmin;							  // index with smallest reliability
	  
	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check
		 double Lnmin2m = DBL_MAX;
		 double Lm2nmin = DBL_MAX;
		 double Lm2star;
		 int alphacum = 0;
		 for(l = 0; l < Nmlen[m]; l++) { // for each adjacent variable node
			n = Nm[m][l];
			alphacum = alphacum ^ alphantom[na[n]][n];
			if(betantom[na[n]][n] < Lnmin2m) { // new smallest incoming
			   Lm2nmin = rcbpboxplus(Lm2nmin,Lnmin2m);
			   nmin = n;
			   Lnmin2m = betantom[na[n]][n];
			}
			else { // not a new min
			   Lm2nmin = rcbpboxplus(Lm2nmin,betantom[na[n]][n]);
			}
		 } // for
		 Lm2star = rcbpboxplus(Lm2nmin,Lnmin2m);

		 // compute the outgoing messages
		 for(l = 0; l < Nmlen[m]; l++) { // for each adjacent variable node
			n = Nm[m][l];
			alphasign = alphacum ^ alphantom[na[n]][n];
			if(n == nmin) { // message to least reliable
			   aminstarLmton[n][na[n]++] = (1 - 2*alphasign)*Lm2nmin;
			}
			else {
			   aminstarLmton[n][na[n]++] = (1 - 2*alphasign)*Lm2star;
			}
		 }
	  } // for m
		 

      if(printstuff) printsparsetranspose("aminstarcorrectLmton",aminstarLmton);

	  // Variable node to check node computation
	  for(n = 0; n < N; n++) {  // work over each column
		 forbacksum(aminstarLc[n],aminstarLmton[n],Mnlen[n],poutsum,sumall);
		 aminstarLcout[n] = sumall;
		 if(aminstarLcout[n] >= 0) c[n] = 0; else c[n] = 1;
		 for(l = 0; l < Mnlen[n]; l++) {
			betantom[l][n] = abs(poutsum[l]);
			alphantom[l][n] = poutsum[l] < 0;
		 }
	  } // for n
	  if(printstuff) printsparse("betantom",betantom);
	  if(printstuff) printsparsechar("alphantom",alphantom);
	  if(printstuff) VECDUMP(MScorrectLcout,N);
	  if(printstuff) VECDUMP2(c,N,int);

	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
		 if(parity) { // if parities check
			numloops = loopcount;
			return parity;
		 }
	  }
   } while(loopcount < maxnumloop);
   if(numloops) { // parity has not been previously checked
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}

void LDPCDECODER::bitflip1decodeinit(double *y, unsigned char *c)
{
   int n;
   int nnonzero = 0;
   for(n = 0; n < N; n++) {
	  c[n] = (y[n] <= 0) ? 0 : 1;
	  if(c[n]) nnonzero++;
   }
}


int LDPCDECODER::bitflip1decode(unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)
{
   int m,n,l;
   unsigned char s1;  // parity sum
   int f1;   // number of unmet paritys sum
   int f1max;
   int flistctr;
   int maxbit;
   int parity = 0;
   bool allparitiesgood;
   unsigned long int loopcount = 0;
   int nnonzerosynd;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
	  allparitiesgood = true;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array
	  int nmin;							  // index with smallest reliability
	  flistctr = 0;						  // number of bits with max parity hit
	  nnonzerosynd = 0;
	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check
		 unsigned char s1 = 0;
		 // compute cH' (mod 2): evaluate parities
		 for(l = 0; l < Nmlen[m]; l++) { // for each bit in the check
			n = Nm[m][l];
			s1 ^= c[n];
		 }
		 synd[m] = s1;
		 if(s1) {
			allparitiesgood = false;  // found a bad syndrome
			nnonzerosynd++;
		 }
	  } // for m
	  if(allparitiesgood) {
		 numloops = loopcount;
		 parity = 1;
		 return parity;  // return 1 = all parities good
	  }
	  // compute s*H (integer): evaluate number of checks affected by each bit
	  f1max = 0;
	  for(n = 0; n < N; n++) {  // for each bit
		 f1 = 0; // number of checks affected by bit n
		 for(l = 0; l < Mnlen[n]; l++) {  // for each check on this bit
			m = Mn[n][l];
			if(synd[m]) {   // syndrome bit set
			   f1++;
			}
		 }
		 if(f1 > f1max) {
			// found a new max --- 
			f1max = f1;
			flistctr = 0;
			flist[flistctr++] = n;
		 }
		 else if(f1 == f1max) {  // another one that has this maximum value
			flist[flistctr++] = n;
		 }
	  } // for n
	  // flip all the bits with the most affected checks
	  for(int i = 0; i < flistctr; i++) {
		 c[flist[i]] = 1-c[flist[i]];
	  }

   } while(loopcount < maxnumloop);

   parity = allparitiesgood;
   numloops = loopcount;
   return parity;
}

void LDPCDECODER::galAbitflipdecodeinit(double *y)
   // convert to posterior probabilities
{
   int i,m,n;
   if(printstuff) VECDUMP(y,N);
   double L;

   // yc = bit version of inputs
   

   // compute channel posterior probabilities
   for(i = 0; i < N; i++) {
	  yc[i] = y[i] < 0 ? 0 : 1;
   }
   // VECDUMP2(yc,N,int);

   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 Bntom[m][n] = yc[n];
	  }
   }
   //    if(printstuff)
   //  printsparsechar("Bntom(init)",Bntom);
}



int LDPCDECODER::galAbitflipdecode(unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)
{
   int m,n,l;
   unsigned char s1;  // parity sum
   int f1;   // number of unmet paritys sum
   int f1max;
  int maxbit;
   int parity = 0;
   bool allparitiesgood;
   unsigned long int loopcount = 0;
   bool allcheckssame;
   unsigned char check;
   int colsum;
   unsigned char bmtoncol;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
	  allparitiesgood = true;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array
	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check

		 unsigned char s1 = 0, s1n;
		 for(l = 0; l < Nmlen[m]; l++) { // for each bit in the check
			n = Nm[m][l];
			s1 ^= Bntom[na[n]][n];
		 }
		 if(s1) allparitiesgood = false;  // found a bad syndrome

		 // for each connected VN, remove its bit influence
		 for(l = 0; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			s1n = s1 ^ Bntom[na[n]][n];
			Bmton[na[n]++][n] = s1n;
		 } // for l
	  } // for m
	  //printsparsechar("Bmton",Bmton);

	  if(allparitiesgood) {
		 numloops = loopcount;
		 return 1;  // return 1 = all parities good
	  }
	  f1max = 0;
	  for(n = 0; n < N; n++) {  // for each bit
		 bntomsum[n] = 0;
		 for(l = 0; l < Mnlen[n]; l++) { // loop over different m indices
			m = Mn[n][l];
			// look at all of the incoming messages (Bmtnon[n][(all - m)])
			allcheckssame = true;
			colsum = 0;  // used to see if all messages are the same
			for(int lp = 0; lp < Mnlen[n]; lp++) {
			   if(lp == l) continue; // check all the others
			   bmtoncol = Bmton[lp][n];   // save in case allcheckssame
			   colsum += Bmton[lp][n];
			} // for l'
			if(!(colsum == (Mnlen[n]-1) || colsum == 0)) { // not all the same
			   allcheckssame = false;
			}
			// if all the incoming messages are the same and
			// the are different from the incoming channel message yc[n],
			// the message out is the complements of yc[n]
			if(allcheckssame && yc[n] != bmtoncol) {
			   Bntom[l][n] = 1-yc[n];
			}
			// else messages don't all agree, or they all agree and equal yc[n]
			else {
			   Bntom[l][n] = yc[n];
			}
			bntomsum[n] += Bntom[l][n]; 
		 }  // for l
		 // Decoding rule:
		 // If the degree of a variable node is odd, the decision is set to
		 // the majority of the incoming check node messages.
		 // If the degree of a variable node is even, the decision is set to
		 // the majority of the majority of (check nodes messages and channel
		 // message)
		 if(Mnlen[n] %2 == 1) {  // degree of the variable node is odd
			if(bntomsum[n] > (Mnlen[n] >> 1))
			   c[n] = 1;
			else
			   c[n] = 0;
		 }
		 else { // the degree of the variable nodes is even
			bntomsum[n] += yc[n];   // include the incoming channel message
			if(bntomsum[n] > ((Mnlen[n]+1) >> 1))
			   c[n] = 1;
			else
			   c[n] = 0;
		 }
	  } // for n
	  
	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
		 if(parity) { // if parities check
			numloops = loopcount;
			return parity;
		 }
	  }
	  if(sigflag) {
		 cout << "sigflag set: breaking" << endl;
		 break;
	  }
   } while(loopcount < maxnumloop);
   if(numloops) { 
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}


void LDPCDECODER::weightedbitflipdecodeinit(double *y)
   // convert to posterior probabilities
{
   int i,m,n;
   if(printstuff) VECDUMP(y,N);

   // VECDUMP(y,N);   
   // compute signed values of received values
   for(i = 0; i < N; i++) {
	  x[i] = y[i] < 0 ? -1 : 1;
x[i] = -x[i];
   }
   // VECDUMP2(x,N,int);
}

int LDPCDECODER::weightedbitflipdecode(double *y, unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)
{
   int m,n,l, l1, np;
   int parity = 0;
   bool allparitiesgood;
   unsigned long int loopcount = 0;
   bool allcheckssame;
   double deltawbfmin;
   double betam;
   double ay;
   int prodx;
   int nmin;

   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
	  allparitiesgood = true;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // compute all the parities and the betas
	  for(m = 0; m < M; m++) {
		 prodx = 1;
		 betam = DBL_MAX;
		 for(l = 0; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			prodx *= x[n];
			ay = abs(y[n]);
			if(ay< betam) betam = ay;
		 }
		 parities[m] = prodx;
		 betas[m] = betam;
		 if(prodx != 1) {
			allparitiesgood = false;
		 }
	  } // for m
      if(allparitiesgood) { // we are done
		 numloops = loopcount;
		 // VECDUMP2(x,N,int);
         for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
		 // VECDUMP2(c,N,int);
         return 1;
	  }

	  // compute the WBF function
	  deltawbfmin = DBL_MAX;
	  for(n = 0; n < N; n++) {
		 deltawbf[n] = 0;
		 for(l = 0; l < Mnlen[n]; l++) {  // for m in M(n)
			m = Mn[n][l];
			deltawbf[n] += betas[m]*parities[m];
		 } // for l (values of m in M(n))
		 if(deltawbf[n] < deltawbfmin) {
			deltawbfmin = deltawbf[n];
			nmin = n;
		 }
	  } // for n
	  if(allparitiesgood) {
		 numloops = loopcount;
		 // convert x to c
		 for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
		 return 1;  // return 1 = all parities good
	  }

	  // flip a bit
	  x[nmin] = -x[nmin];
// VECDUMP2(x,N,int);
	  if(sigflag) {
		 cout << "sigflag set: breaking" << endl;
		 break;
	  }
   } while(loopcount < maxnumloop);
   parity = 0; // if we get to here, parities didn't check
   numloops = loopcount;
   for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
   return parity;
}





void LDPCDECODER::modweightedbitflipdecodeinit(double *y)
   // convert to posterior probabilities
{
   int i,m,n;
   if(printstuff) VECDUMP(y,N);

   // VECDUMP(y,N);   
   // compute signed values of received values
   for(i = 0; i < N; i++) {
	  x[i] = y[i] < 0 ? -1 : 1;
x[i] = -x[i];
   }
   // VECDUMP2(x,N,int);
}

int LDPCDECODER::modwtedbitflipdecode(double *y, unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)
{
   int m,n,l, l1, np;
   int parity = 0;
   bool allparitiesgood;
   unsigned long int loopcount = 0;
   bool allcheckssame;
   double deltawbfmin;
   double betam;
   double ay;
   int prodx;
   int nmin;
   int sumparities;

   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
	  allparitiesgood = true;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // compute all the parities and the betas
	  sumparities = 0;
	  for(m = 0; m < M; m++) {
		 prodx = 1;
		 betam = DBL_MAX;
		 for(l = 0; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			prodx *= x[n];
			ay = abs(y[n]);
			if(ay< betam) betam = ay;
		 }
		 parities[m] = prodx;
		 sumparities += prodx;
		 betas[m] = betam;
		 if(prodx != 1) {
			allparitiesgood = false;
		 }
	  } // for m
// printf("wbf: loopcount=%ld  sumparities=%d\n",loopcount,sumparities);
      if(allparitiesgood) { // we are done
		 numloops = loopcount;
		 // VECDUMP2(x,N,int);
         for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
		 // VECDUMP2(c,N,int);
// printf("wbf: loopcount=%ld  return=0\n",loopcount);
         return 1;
	  }

	  // compute the WBF function
	  deltawbfmin = DBL_MAX;
	  for(n = 0; n < N; n++) {
		 deltawbf[n] = modweightedbitflipalpha*abs(y[n]);
		 for(l = 0; l < Mnlen[n]; l++) {  // for m in M(n)
			m = Mn[n][l];
			deltawbf[n] += betas[m]*parities[m];
		 } // for l (values of m in M(n))
		 if(deltawbf[n] < deltawbfmin) {
			deltawbfmin = deltawbf[n];
			nmin = n;
		 }
	  } // for n
	  if(allparitiesgood) {
		 numloops = loopcount;
		 // convert x to c
		 for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
		 return 1;  // return 1 = all parities good
	  }

	  // flip a bit
	  x[nmin] = -x[nmin];
// VECDUMP2(x,N,int);
	  if(sigflag) {
		 cout << "sigflag set: breaking" << endl;
		 break;
	  }
   } while(loopcount < maxnumloop);
   parity = 0; // if we get to here, parities didn't check
// printf("wbf: loopcount=%ld  return=0\n",loopcount);
   numloops = loopcount;
   for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
   return parity;
}
 

void LDPCDECODER::grad1bitflipdecodeinit(double *y)
   // convert to posterior probabilities
{
   int i,m,n;
   if(printstuff) VECDUMP(y,N);

   // VECDUMP(y,N);   
   // compute signed values of received values
   for(i = 0; i < N; i++) {
	  x[i] = y[i] < 0 ? -1 : 1;
x[i] = -x[i];
   }
   // VECDUMP2(x,N,int);
}

int LDPCDECODER::grad1bitflipdecode(double *y, unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)
{
   int m,n,l, l1, np;
   int parity = 0;
   bool allparitiesgood;
   unsigned long int loopcount = 0;
   bool allcheckssame;
   double deltawbfmax;
   double betam;
   double ay;
   char prodx;
   int nmax;

   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   // VECDUMP2(x,N,int);
   do {
	  loopcount++;
	  allparitiesgood = true;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // compute all the parities and the betas
	  for(m = 0; m < M; m++) {
		 prodx = 1;
		 for(l = 0; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			prodx *= x[n];
		 }
		 parities[m] = prodx;
		 if(prodx != 1) {
			allparitiesgood = false;
		 }
	  } // for m
	  // VECDUMP2(parities,M,int);

      if(allparitiesgood) { // we are done
		 numloops = loopcount;
		 // VECDUMP2(x,N,int);
         for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
		 // VECDUMP2(c,N,int);
         return 1;
	  }

	  // compute the GD function
	  deltawbfmax = DBL_MIN;
	  for(n = 0; n < N; n++) {
		 deltawbf[n] = x[n]*y[n];
		 for(l = 0; l < Mnlen[n]; l++) {  // for m in M(n)
			m = Mn[n][l];
			deltawbf[n] += parities[m];
		 } // for l (values of m in M(n))
		 if(abs(deltawbf[n]) > deltawbfmax) {
			deltawbfmax = abs(deltawbf[n]);
			nmax = n;
		 }
	  } // for n
	  // VECDUMP(deltawbf,N);
	  // cout << "deltawbfmax=" << deltawbfmax << "  nmax=" << nmax << endl;
	  // flip a bit
	  x[nmax] = -x[nmax];
	  // VECDUMP2(x,N,int);
	  if(sigflag) {
		 cout << "sigflag set: breaking" << endl;
		 break;
	  }
   } while(loopcount < maxnumloop);
   parity = 0; // if we get to here, parities didn't check
   numloops = loopcount;
   for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
   return parity;
}

// grad2bitflipdecode:  multi-bit gradient descent bit flip

void LDPCDECODER::grad2bitflipdecodeinit(double *y)
   // convert to posterior probabilities
{
   int i,m,n;
   if(printstuff) VECDUMP(y,N);

   // compute signed values of received values
   for(i = 0; i < N; i++) {
	  x[i] = a_sign*y[i] >= 0 ? -1 : 1;
	  //  x[i] = -x[i];
   }

   // set decoder parameters
   grad2bitflip_mu = 0;			// start decoding with multibit flipping
   grad2bitflip_theta = -0.6;   // flipping threshold

   fmin = DBL_MAX;    // find range of function values
   fmax = -DBL_MAX;
// cout << "deltawbfmin=" << deltawbfmin << "    deltawbfmax=" << deltawbfmax << endl;
}

int LDPCDECODER::grad2bitflipdecode(double *y, unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)
{
   int m,n,l, l1, np;
   int parity = 0;
   bool allparitiesgood;
   unsigned long int loopcount = 0;
   bool allcheckssame;
   double betam;
   double ay;
   char prodx;
   double f1, f2;
   int sumprodx;
   int numtoflip;
   int nminidx;
   int sumparities1, sumparities2;
   int maxsumparities = 0;
   bool docheckparity;

   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   docheckparity = true;
   // VECDUMP2(x,N,int);
   do {
	  deltawbfmin = DBL_MAX;
	  deltawbfmax = -DBL_MAX;

	  loopcount++;
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  if(docheckparity) {
		 gdcomputeparities(x,parities,sumparities1,allparitiesgood);
		 docheckparity = false;
	  }

      if(allparitiesgood) { // we are done
		 numloops = loopcount;
         for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
// printf("exiting: loopcount=%ld  return=1\n",loopcount);
         return 1;
	  }
	  // compute the objective function and the GD function

	  f1 = gdobjectivefunct(y,x, parities, sumparities1, deltawbf,
							deltawbfmin, deltawbfmax, nminidx);
// printf("mbgd: mu=%d  loopcount=%ld  sumparities1=%d  f1=%g\n",int(grad2bitflip_mu), loopcount,sumparities1,
// 			 f1);

	  numtoflip = 0;

	  if(!grad2bitflip_mu) {   // multibit mode
		 // find the smallest grad2bitflip_percent of the deltawbf values
		 // and flip the bits there
		 double grad2bitflip_percent = .05;
		 double deltarange, deltathresh;
		 deltarange = deltawbfmax - deltawbfmin;
		 deltathresh = deltawbfmin + grad2bitflip_percent*deltarange;
		 if(deltathresh > 0) deltathresh = -.1;
		 numtoflip = 0;


		 // try flipping all of these bits
		 for(n = 0; n < N; n++) {
			if(deltawbf[n] <= deltathresh) {
			   numtoflip++;
			   tryx[n] = 1-x[n];
			}
			else tryx[n] = x[n];
		 }

		 // compute the objective function again
		 gdcomputeparities(tryx,parities,sumparities2,allparitiesgood);
		 f2 = gdobjectivefunct(y,tryx, parities, sumparities2, deltawbf,
							   deltawbfmin, deltawbfmax, nminidx);

// printf("mbgd: mu=%d  f1=%g  f2=%g  sumparities1=%d  sumparities2=%d\n",
// int(grad2bitflip_mu),f1,f2,sumparities1, sumparities2);
		 if(f2 < f1) {					   // things got worse
			grad2bitflip_mu = 1; // switch to single flip mode
			// don't keep the flip
// printf("### switching to single flip mode\n");
            docheckparity = true;
		 }
		 else { // things got better -- keep the current "try" bit flips
			char *tempx;
			tempx = x;
			x = tryx;
			tryx = tempx;
			docheckparity = false;  // already checked
			sumparities1 = sumparities2;
		 }
	  }
	  else {  // single bit
// printf("singlebitflip: minidx=%d   deltawbfmin=%g\n",nminidx,deltawbfmin);
		 x[nminidx] = -x[nminidx];

		 docheckparity = true;
// 		 gdcomputeparities(x,parities,sumparities2,allparitiesgood);
// 		 docheckparity = false;

// docheckparity = true;

// // for kicks
// f2= gdobjectivefunct(y,x, parities, sumparities2, deltawbf,
// 							deltawbfmin, deltawbfmax, nminidx);


// printf("mu=%d  deltawbfmin=%g  sumparities2=%d   f2=%g\n",grad2bitflip_mu, deltawbfmin,sumparities2,f2);

//           sumparities1 = sumparities2;
	  }
		 
	  // Gather some statistics on f1, f2  (for debug/test)
	  if(f1 < fmin){
		 fmin = f1;
	  }
	  if(f1 > fmax) {
		 fmax = f1;
	  }
	  
	  // VECDUMP2(x,N,int);
	  if(sigflag) {
		 cout << "sigflag set: breaking" << endl;
		 break;
	  }
   } while(loopcount < maxnumloop);

   parity = 0; // if we get to here, parities didn't check
// printf("mbgd: loopcount=%ld  return=0\n",loopcount);
   numloops = loopcount;
   for(n = 0; n < N; n++) c[n] = (1 - x[n])/2;
// printf("exiting: parity=%d  f1=%g  f2=%g  fmin=%g  fmax=%g\n",int(parity),
// f1,f2,fmin,fmax);
   return parity;
}


void LDPCDECODER::gdcomputeparities(char *x, char *parities,
									int &sumparities,
									  bool &allparitiesgood)
{
   allparitiesgood = true;
   int prodx = 1;
   sumparities = 0;
   int n;
   
   for(int m = 0; m < M; m++) {
	  prodx = 1;
	  for(int l = 0; l < Nmlen[m]; l++) {
		 n = Nm[m][l];
		 prodx *= x[n];
	  }
	  parities[m] = prodx;
	  sumparities += prodx;
	  if(prodx != 1) {
		 allparitiesgood = false;
	  }
   } // for m
}


double LDPCDECODER::gdobjectivefunct(double *y,char *x, char *parities,
									 int &sumparities,
									 double *deltawbf,
									 double &deltawbfmin,
									 double &deltawbfmax,
									 int &nminidx)
// x here overshadows the x defined in the class
// sumparities and parities are input variables;
// they must be computed prior to calling this function
// 
{
   int m, l, n;
   char prodx;
   double f;
   double xy;
   
   f = sumparities;
   for(n = 0; n < N; n++) {
	  xy = x[n]*(a_sign*y[n]);
	  f += xy;
	  for(l = 0; l < Mnlen[n]; l++) {  // for m in M(n)
		 m = Mn[n][l];
		 deltawbf[n] += parities[m];
	  } // for l (values of m in M(n))
	  if(deltawbf[n] < deltawbfmin) {
		 deltawbfmin = deltawbf[n];
		 nminidx = n;
	  }
	  if(deltawbf[n] > deltawbfmax) {
		 deltawbfmax = deltawbf[n];
	  }
	  
   }
   return f;
}
	  
// --------------------------------------------------------

void LDPCDECODER::DC1decodeinit(double *y)
{
   int n,m;
   double Lc1;
   double Labssum = 0;
   double DC1eps = 0.01;   		// settable parameter
   double xp,p;

   Lc1 = 2*a/sigma2;// a is negative if the bit 0 is represented with a negative

   Lcsumsq = 0;
   for(n = 0; n < N; n++) {
	  Lc[n] = Lc1*y[n];   // Lc[n]>0: bit probably 0
	  Lcsumsq += Lc[n]*Lc[n];
	  Labssum += abs(Lc[n]);
   }
if(debugprint){ VECDUMP(y,N); VECDUMP(Lc,N); }
 
   // set the Emax constraint
   EmaxDC = -(1.0 + DC1eps)*Labssum;
   if(debugprint) { cout << "EmaxDC=" << EmaxDC << endl;
	  cout << "DCdoprobconstraint=" << int(DCdoprobconstraint) << endl; }

   if(DCdoprobconstraint) {  // initialize with: p = exp(Lc[n])/(1+exp(Lc[n]),
	  for(n = 0; n < N; n++) {
		 xp = exp(Lc[n])/(1 + exp(Lc[n]));   // Prob(X=1)
		 p = 2*xp - 1;   // p>0: X probably 1.  p<0: X probably 0.
if(debugprint) {cout << "Lc=" << Lc[n] << "  xp=" << xp << "  2xp-1=" << p << "  ";}
		 for(m = 0; m < Mnlen[n]; m++) {
			DCrdat1[m][n] = p;
		 }
		 DCrdat1[maxcolwt][n] = p;
	  }
if(debugprint) {cout << endl;}
   }
   else {   // use likelihood data to initialize
	  for(n = 0; n < N; n++) {
		 for(m = 0; m < Mnlen[n]; m++) {
			DCrdat1[m][n] = Lc[n];
		 }
		 DCrdat1[maxcolwt][n] = Lc[n];
	  }
   }
// cout << "a=" << a << "  a_sign=" << a_sign << "  Lc1=" << Lc1 << endl;
if(debugprint) {printsparseMp1("DC1",DCrdat1);}
// cout << "end initialize" << endl;
}
   

int LDPCDECODER::DC1decode(double * y, unsigned char *c,
							 unsigned long int maxnumloop,
							 unsigned long int &numloops)
{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;
   double **temppt;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }
   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array
	  parity = DC1decodePD(DCrdat1,DCrdat2);
if(debugprint) {
cout << "After PD:" << endl;
printsparseMp1("DCrdat1:",DCrdat1);
 printsparseMp1("PD (DCrdat2):",DCrdat2); }

	  DC1DM1(DCrdat1, DCrdat2);
	  // DCrdat2 <- DCrdat2 - DCrat1
	  // DCrat1 <- 2*DCrat2 - DCrat1
if(debugprint) printsparseMp1("Pd(rt) - rt",DCrdat2);
if(debugprint) printsparseMp1("2 Pd(rt) - rt",DCrdat1);
	  DC1decodePC(DCrdat1, DCrdat1, c);
	  // DCrdat1 = Pc(2PD(rt) - rt);
if(debugprint) printsparseMp1("Pc(2 Pd(rt) - rt)",DCrdat1);
	  parity = checkparity(c);
if(debugprint) { cout << "parity=" << parity << endl;}
	  if(parity) {
// cout << "done!" << endl;
		 numloops = loopcount;
cout << "numloops=" << numloops << endl;
		 return parity;
	  }
	  DC1DM2(DCrdat1,DCrdat2);
	  // DCrdat1 = DCrdat1 - DCrdat2
if(debugprint) printsparseMp1("Next set of parameters",DCrdat1);
   } while(loopcount < maxnumloop);
   numloops = loopcount;
   return parity;
}

int LDPCDECODER::DC1decodePD(double **rin, double **rout)
// Do the Divide projection on the input data (rin)
// to produce the output data (rout)
// This cannot be done in place (rout cannot equal rin)
{
   int m, l, n, nuidx, nminus;
   double numin, ah, rmn;
   bool allparitiesgood;
   double Lx;
   int parity;
   
   bzero(na,N*sizeof(unsigned int));

   // Divide Step
   allparitiesgood = true;
   // project onto the parity constraints
   for(m = 0; m < M; m++) {
	  nminus = 0;
	  numin = DBL_MAX;
	  for(l = 0; l < Nmlen[m]; l++) {  // for each replica on this row
		 n = Nm[m][l];
		 rmn = rin[na[n]][n];
		 x[n] = (rmn >= 0 ? 1 : -1);  // x[n]=1 corresponds to bit=0
		 ah = abs(rmn);
		 if(ah < numin) {
			numin = ah;
			nuidx = n;
		 }
		 if(x[n] < 0) nminus++;
	  }
	  if(nminus % 2) {		// odd number of -1s --- flip the smallest
if(debugprint) {cout << "flipping bit: m=" << m << "  nuidx=" << nuidx << "  numin=" << numin << endl;}
		 x[nuidx] = -x[nuidx];
		 allparitiesgood = false;
	  }
	  // copy the projected data into the output array
	  for(l = 0; l < Nmlen[m]; l++) {
		 n = Nm[m][l];
		 rout[na[n]++][n] = x[n];
	  }
   } // for m

   if(allparitiesgood) {
	  parity = 1;
	  return parity;  // return 1 = all parities good
   }

   // project onto the energy constraint
   m = maxcolwt; // constraint number M
if(debugprint){
cout << "r before energy projection: " << endl;
VECDUMP(rin[m],N);}
   Lx = 0;  // EmaxDC;
   for(n = 0; n < N; n++) {
	  Lx -= Lc[n]*rin[m][n];
   }
if(debugprint) {cout << "energy before projection: " << Lx << endl;}
   if(!(Lx <= EmaxDC)) {
if(debugprint){ cout << "constraint not satisfied" << endl;}
	  Lx = -Lx + EmaxDC;
	  Lx = Lx/Lcsumsq;
	  for(n = 0; n < N; n++) {
		 rout[m][n] = rin[m][n] - Lc[n]*Lx;
	  }
if(debugprint) {cout << "r after energy projection: " << endl;
VECDUMP(rout[m],N);}
   }
   else {
if(debugprint) {cout << "energy constraint already satisfied" << endl;}
      m = maxcolwt; // constraint number M
	  for(n = 0; n < N; n++) {
		 rout[m][n] = rin[m][n] - Lc[n]*Lx;
	  }
   }
// Check satisfaction of constraint   
// double xlsum = 0;
//  cout << "x: ";
// VECDUMP2(x,N,int);
// for(n = 0; n < N; n++) {
//    xlsum += Lc[n]*x[n];
// }
// cout << "xlsum = " << xlsum << "  EmaxDC=" << EmaxDC << endl;

// Lx = -EmaxDC;
// for(n = 0; n < N; n++) {
// Lx -= Lc[n]*rout[m][n];
// }
// cout << "constraint value=" << Lx << endl;
   // (should be < 0)

   parity = 0;
   return parity;
}

void LDPCDECODER::DC1decodePC(double **rin, double **rout, unsigned char *c)
// Do the Concur projection on the input data (rin)
// to produce the output data (rout)
//
// This can be done in place (rin can equal rout)
{
   int n,m;
   double rave;
   // Concur Step

if(debugprint){printsparseMp1("PCin",rin);}

   for(n = 0; n < N; n++) {
	  rave = rin[maxcolwt][n];   // the energy constraint
	  for(m = 0; m < Mnlen[n]; m++) {
		 rave += rin[m][n];
	  }
	  rave /= (Mnlen[n] + 1);
	  for(m = 0; m < Mnlen[n]; m++) {
		 rout[m][n] = rave;
	  }
	  rout[maxcolwt][n] = rave;
	  if(rave > 0) {
		 c[n] = 0;
	  }
	  else {
		 c[n] = 1;
	  }
   }
if(debugprint) {printsparseMp1("PCout",rout);
VECDUMP2(c,N,int);}
}

void LDPCDECODER::DC1DM1(double **DCrdat1, double **DCrdat2)
// Input:
//  DCrdat1 = r_t
//  DCrdat2 = PD(r_t)
//
// Output
// DCrdat2 = PD(r_t) - r_t
// DCrdat1 = 2*PD(r_t) - r_t
{

if(debugprint){
cout << "DC1DM1" << endl;
printsparseMp1("DCrdat1",DCrdat1);
printsparseMp1("DCrdat2",DCrdat2);}

   int n,m;
   double t;
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 t = DCrdat2[m][n];  // Pd(rt)
		 DCrdat2[m][n] = t - DCrdat1[m][n];
		 DCrdat1[m][n] = 2*t - DCrdat1[m][n];
	  }
	  t = DCrdat2[maxcolwt][n];
	  DCrdat2[maxcolwt][n] = t - DCrdat1[maxcolwt][n];
	  DCrdat1[maxcolwt][n] = 2*t - DCrdat1[maxcolwt][n];
   }
// if(debugprint) {
// printsparseMp1("Pd(r) - r",DCrdat2);
// printsparseMp1("2*Pd(r) - r",DCrdat1);}

}

void LDPCDECODER::DC1DM2(double **DCrdat1, double **DCrdat2)
// Input:
// DCrdat1 = Pc(2*PD(r_t) - r_t)
// DCrdat2 = PD(r_t) - r_t
//
// Output:   
// DCrdat1 = Pc(2*Pd(r_t) - r_t - (Pd(r_t) - r_t)
// DCrdat2 = PD(r_t) - r_t
{
   int n,m;
   
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 DCrdat1[m][n] = DCrdat1[m][n] - DCrdat2[m][n];
	  }
	  m = maxcolwt;
	  DCrdat1[m][n] = DCrdat1[m][n] - DCrdat2[m][n];
   }
// if(debugprint){printsparseMp1("Pc(2*Pd(r) - r) - (Pd(r_t) - r_t)",DCrdat1);}
}

// -------------------------------------------------------------------

void LDPCDECODER::DMBPdecodeinit(double *y)
{
   int n,m;
   double Lc;
   Lc = 2*a/sigma2;
   for(n = 0; n < N; n++) {
	  MSLc[n] = Lc*y[n];
   }

   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 alphantom[m][n] = (MSLc[n] < 0);

		 betantom[m][n] = abs(MSLc[n]);
	  }
   }
   if(printstuff) printsparse("DMBP correct Lntom(beta,init)",betantom);
   if(printstuff) printsparsechar("alphantom",alphantom);
}

int LDPCDECODER::DMBPdecode(unsigned char *c,
									 unsigned long int maxnumloop,
									 unsigned long int &numloops)

{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;
   char alphasign, alphacum;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check

		 // work over nonzero elements of this row
		 l = 0;
		 n = Nm[m][l];
		 alphasign = alphantom[na[n]][n];
		 alpharow[l] = alphasign;
		 alphacum = alphasign;
		 betarow[l] = betantom[na[n]][n];
		 for(l = 1; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			alphasign = alphantom[na[n]][n];
			alphacum ^= alphasign;
			alpharow[l] = alphasign;
			betarow[l] = betantom[na[n]][n];
		 }
		 // compute leave-one-out min data
		 forbackmincorrect(betarow,Nmlen[m],poutmin);
		 for(l = 0; l < Nmlen[m]; l++) { // for each n on this row
			n = Nm[m][l];
			alphasign = alphacum ^ alphantom[na[n]][n]; 
			MScorrectLmton[n][na[n]++] = (1-2*alphasign)*poutmin[l];
		 }
	  }  // for m
      if(printstuff) printsparsetranspose("MScorrectLmton",MScorrectLmton);
	  if(savestuff1) savestuff(MScorrectLmton,fsave2);

	  // Update beliefs (at each variable node)
	  for(n = 0; n < N; n++) {  // work over each column
		 DMBPb[n] = MSLc[n];   // channel input
		 for(l = 0; l < Mnlen[n]; l++) {  // no leave-one-out here
			DMBPb[n] += MScorrectLmton[n][l];
		 }
		 DMBPb[n] *= DMBPZ;
		 c[n] = (DMBPb[n] < 0? 1: 0);
	  }
	  parity = checkparity(c);
	  if(parity) {
		 numloops = loopcount;
		 return parity;
	  }

	  // Update messages from variables to checks
	  double Mntom;
	  for(n = 0; n < N; n++) {
		 for(l = 0; l < Mnlen[n]; l++) {
			Mntom = DMBPb[n] - 0.5*(MScorrectLmton[n][l] -
									(1-2*alphantom[l][n])*betantom[l][n]);
			alphantom[l][n] = (Mntom < 0);
			betantom[l][n] = abs(Mntom);
		 }
	  }
	  
   } while(loopcount < maxnumloop);
   numloops = loopcount;
   return parity;
}


// --------------------------------------------------------------
// Utility functions

void LDPCDECODER::printsparse(string name,double **sparsemat)
{
   int l,n, m,lastn,n1;
   int *na;

   int sp = cout.precision();
   cout.precision(2);
   na = new int[N];
   for(l = 0; l < N; l++) na[l] = 0;
   cout << name << endl;
   for(m = 0; m < M; m++) {
	  lastn = 0;
	  for(l = 0; l < Nmlen[m]; l++) {
		 n = Nm[m][l];
		 // print 0s as necessary
		 for(n1 = lastn; n1 < n; n1++) {
			cout << setw(5) << "0\t";
		 }
		 // print the new number
		 cout << setw(5) << sparsemat[na[n]++][n] << "\t";
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < N; n1++) {
		 cout << setw(5) << "0.0\t";
	  }
	  cout << endl;
   }
   delete[] na;
   cout.precision(sp);
}

void LDPCDECODER::printsparseMp1(string name,double **sparsemat)
// print a sparse matrix, and one more row (to show energy constraints)
{
   int n;

   int sp = cout.precision();
   cout.precision(2);
   printsparse(name,sparsemat);
   cout << endl;  // skip a line
   for(n = 0; n < N; n++) {
	  cout << setw(5) << sparsemat[maxcolwt][n] << "\t";
   }
   cout.precision(sp);
   cout << endl;
}
   
void LDPCDECODER::printsparsechar(string name,char **sparsemat)
{
   int l,n, m,lastn,n1;
   int *na;

   int sp = cout.precision();
   cout.precision(2);
   na = new int[N];
   for(l = 0; l < N; l++) na[l] = 0;
   cout << name << endl;
   for(m = 0; m < M; m++) {
	  lastn = 0;
	  for(l = 0; l < Nmlen[m]; l++) {
		 n = Nm[m][l];
		 // print 0s as necessary
		 for(n1 = lastn; n1 < n; n1++) {
			cout << setw(5) << "0";
			cout << "\t";
		 }
		 // print the new number
		 cout << setw(5) << int(sparsemat[na[n]++][n]);
		 cout << "\t";
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < N; n1++) {
		 cout << setw(5) << "0\t";
	  }
	  cout << endl;
   }
   delete[] na;
   cout.precision(sp);
}

void LDPCDECODER::printsparsechar(string name,unsigned char **sparsemat)
{
   int l,n, m,lastn,n1;
   int *na;

   int sp = cout.precision();
   cout.precision(2);
   na = new int[N];
   for(l = 0; l < N; l++) na[l] = 0;
   cout << name << endl;
   for(m = 0; m < M; m++) {
	  lastn = 0;
	  for(l = 0; l < Nmlen[m]; l++) {
		 n = Nm[m][l];
		 // print 0s as necessary
		 for(n1 = lastn; n1 < n; n1++) {
			cout << setw(5) << "0\t";
		 }
		 // print the new number
		 cout << setw(5) << int(sparsemat[na[n]++][n]) << "\t";
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < N; n1++) {
		 cout << setw(5) << "0\t";
	  }
	  cout << endl;
   }
   delete[] na;
   cout.precision(sp);
}


void LDPCDECODER::printsparsetranspose(string name,double **sparsemat)
{
   int l,n, m,lastn,n1;
   int *na;

   int sp = cout.precision();
   cout.precision(2);
   na = new int[N];
   bzero(na,N*sizeof(int));
   cout << name << endl;
   for(m = 0; m < M; m++) {
	  lastn = 0;
	  for(l = 0; l < Nmlen[m]; l++) {
		 n = Nm[m][l];
		 // print 0s as necessary
		 for(n1 = lastn; n1 < n; n1++) {
			cout << setw(5) << "0\t";
		 }
		 // print the new number
		 cout << setw(5) << sparsemat[n][na[n]++]<< "\t";
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < N; n1++) {
		 cout << setw(5) << "0.0\t";
	  }
	  cout << endl;
   }
   delete[] na;
   cout.precision(sp);
}

void
LDPCDECODER::printsparseA(void)  // print the sparse binary A matrix
{
   int l,n, m,lastn,n1;
   int lastm,m1;

   for(m = 0; m < M; m++) {
	  lastn = 0;
	  for(l = 0; l < Nmlen[m]; l++) {
		 n = Nm[m][l];
		 // print 0s as necessary
		 for(n1 = lastn; n1 < n; n1++) {
			cout << "0 ";
		 }
		 // print the new number
		 cout << "1 "; 
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < N; n1++) {
		 cout << "0 ";
	  }
	  cout << endl;
   }
}


inline
void LDPCDECODER::forbackprobs(double *p, int m,double *poutprods)
{
   if(m == 1) {
	  poutprods[0] = 1;
   }
   else {
	  // build the forward/backward tables
	  forwardprod[0] = p[0];           backwardprod[0] = p[m-1];
	  for(int i = 1; i <= m-2; i++) {
		 forwardprod[i] = p[i] * forwardprod[i-1];
		 backwardprod[i] = p[m-1-i]*backwardprod[i-1];
	  }
	  // pull out the partial products
	  poutprods[0] = backwardprod[m-2];
	  for(int i = 1; i < m-1; i++) {
		 poutprods[i] = forwardprod[i-1]*backwardprod[m-2-i];
	  }
	  poutprods[m-1] = forwardprod[m-2];
   }
}

void LDPCDECODER::forbackprobs2(double b0, double *p0, double b1,double *p1,
								int m,
								double *poutprods0,double &pall0,
								double *poutprods1,double &pall1)
{
   // build the forward/backward tables
   if(m == 1) {
	  poutprods0[0] = b0;
	  poutprods1[0] = b1;
	  pall0 = b0*p0[0];
	  pall1 = b1*p1[0];
   }
   else {
	  forwardprod0[0] = b0*p0[0];           backwardprod0[0] = p0[m-1];
	  forwardprod1[0] = b1*p1[0];           backwardprod1[0] = p1[m-1];
	  for(int i = 1; i <= m-2; i++) {
		 forwardprod0[i] = p0[i] * forwardprod0[i-1];
		 backwardprod0[i] = p0[m-1-i]*backwardprod0[i-1];
		 forwardprod1[i] = p1[i] * forwardprod1[i-1];
		 backwardprod1[i] = p1[m-1-i]*backwardprod1[i-1];
	  }
	  // pull out the partial products
	  poutprods0[0] = b0*backwardprod0[m-2];
	  poutprods1[0] = b1*backwardprod1[m-2];
	  for(int i = 1; i < m-1; i++) {
		 poutprods0[i] = forwardprod0[i-1]*backwardprod0[m-2-i];
		 poutprods1[i] = forwardprod1[i-1]*backwardprod1[m-2-i];
	  }
	  poutprods0[m-1] = forwardprod0[m-2];
	  pall0 = forwardprod0[m-2]*backwardprod0[0];
	  poutprods1[m-1] = forwardprod1[m-2];
	  pall1 = forwardprod1[m-2]*backwardprod1[0];
   }
}

void LDPCDECODER::forbacksum(double b0, double *p0, int m,
							 double *poutsum0,double &psum0)
{
   // build the forward/backward tables
   if(m == 1) {
	  poutsum0[0] = b0;
	  psum0 = b0 + p0[0];
   }
   else {
	  forwardsum0[0] = b0 + p0[0];           backwardsum0[0] = p0[m-1];
	  for(int i = 1; i <= m-2; i++) {
		 forwardsum0[i] = p0[i] +  forwardsum0[i-1];
		 backwardsum0[i] = p0[m-1-i] + backwardsum0[i-1];
	  }
	  // pull out the partial products
	  poutsum0[0] = b0 + backwardsum0[m-2];
	  for(int i = 1; i < m-1; i++) {
		 poutsum0[i] = forwardsum0[i-1] + backwardsum0[m-2-i];
	  }
	  poutsum0[m-1] = forwardsum0[m-2];
	  psum0 = forwardsum0[m-2] + backwardsum0[0];
   }
}


void LDPCDECODER::forbackmin(double *p0, int m,double *poutmin0)
{
   // build the forward/backward tables
   if(m == 1) {
	  poutmin0[0] = 0;
   }
   else {
	  forwardmin0[0] = p0[0];           backwardmin0[0] = p0[m-1];
	  for(int i = 1; i <= m-2; i++) {
		 forwardmin0[i] = std::min(p0[i],forwardmin0[i-1]);
		 backwardmin0[i] = std::min(p0[m-1-i],backwardmin0[i-1]);
	  }
	  // pull out the partial products
	  poutmin0[0] = backwardmin0[m-2];
	  for(int i = 1; i < m-1; i++) {
		 poutmin0[i] = std::min(forwardmin0[i-1],backwardmin0[m-2-i]);
	  }
	  poutmin0[m-1] = forwardmin0[m-2];
   }
}

void LDPCDECODER::forbackmincorrect(double *p0, int m,double *poutmin0)
{
   // build the forward/backward tables
   if(m == 1) {
	  poutmin0[0] = 0;
   }
   else {
	  forwardmin0[0] = p0[0];
	  backwardmin0[0] = p0[m-1];
	  for(int i = 1; i <= m-2; i++) {
		 forwardmin0[i] = std::min(p0[i],forwardmin0[i-1]) +
			stildecorrect(p0[i],forwardmin0[i-1]);
		 backwardmin0[i] = std::min(p0[m-1-i],backwardmin0[i-1])
			+ stildecorrect(p0[m-1-i],backwardmin0[i-1]);
	  }
	  // pull out the partial products
	  poutmin0[0] = backwardmin0[m-2];
	  for(int i = 1; i < m-1; i++) {
		 poutmin0[i] = std::min(forwardmin0[i-1],backwardmin0[m-2-i])
			+ stildecorrect(forwardmin0[i-1],backwardmin0[m-2-i]);
	  }
	  poutmin0[m-1] = forwardmin0[m-2];
   }
}

double LDPCDECODER::scorrect(double x,double y)
{
   return log(1 + exp(-fabs(x+y))) - log(1 + exp(-fabs(x - y)));
}

double LDPCDECODER::stildecorrect(double x,double y)
{
   const double c = 0.5;
   double xpy = fabs(x+y);
   double xmy = fabs(x-y);
   if(xpy < 2 && (xmy > 2*xpy)) return c;
   if(xmy < 2 && (xpy > 2*xmy)) return -c;
   return 0;
}


// ------------------------------------------------------------
// Stuff for phi-based decoding
void LDPCDECODER::phidecodeinit(double *y)
{
   int n,m;
   double Lc;
   Lc = 2*a/sigma2;
   for(n = 0; n < N; n++) {
	  phiLc[n] = Lc*y[n];
   }
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 alphantom[m][n] = (phiLc[n] < 0);
		 betantom[m][n] = abs(phiLc[n]);
	  }
   }
   if(printstuff) printsparse("phiLntom(beta,init)",betantom);
   if(printstuff) printsparsechar("alphantom",alphantom);

}


int LDPCDECODER::phidecode(unsigned char *c, unsigned long int maxnumloop,
							unsigned long int &numloops)
{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;
   char alphasign, alphacum;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check
		 // work over nonzero elements of this row
		 l = 0;
		 n = Nm[m][l];
		 alphasign = alphantom[na[n]][n];
		 alpharow[l] = alphasign;
		 alphacum = alphasign;
		 betarow[l] = phi(betantom[na[n]][n]);
		 for(l = 1; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			alphasign = alphantom[na[n]][n];
			alphacum ^= alphasign;
			alpharow[l] = alphasign;
			phibetarow[l] = phi(betantom[na[n]][n]);
		 }
		 // compute leave-one-out min data
		 forbacksum(0,phibetarow,Nmlen[m],phioutphisum,sumall);
		 for(l = 0; l < Nmlen[m]; l++) { // for each n on this row
			n = Nm[m][l];
			alphasign = alphacum ^ alphantom[na[n]][n]; 
			phiLmton[n][na[n]++] = (1-2*alphasign)*phi(phioutphisum[l]);
		 }
	  }  // for m
	  
      if(printstuff) printsparsetranspose("phiLmton",phiLmton);

	  // Variable node to check node computation
	  for(n = 0; n < N; n++) {  // work over each column
		 forbacksum(phiLc[n],phiLmton[n],Mnlen[n],poutsum,sumall);
		 phiLcout[n] = sumall;
		 if(phiLcout[n] >= 0) c[n] = 0; else c[n] = 1;
		 for(l = 0; l < Mnlen[n]; l++) {
			betantom[l][n] = abs(poutsum[l]);
			alphantom[l][n] = poutsum[l] < 0;
		 }
	  } // for n
	  if(printstuff) printsparse("betantom",betantom);
	  if(printstuff) printsparsechar("alphantom",alphantom);
	  if(printstuff) VECDUMP(phiLcout,N);
	  if(printstuff) VECDUMP2(c,N,int);

	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
		 if(parity) { // if parities check
			numloops = loopcount;
			return parity;
		 }
	  }
   } while(loopcount < maxnumloop);
   if(numloops) { // parity has not been previously checked
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}

int LDPCDECODER::qphidecode(unsigned char *c, unsigned long int maxnumloop,
							unsigned long int &numloops)
{
   int m,n1,n2;
   int l,k,idx,n,i;
   int parity = 0;
   double sumall;
   unsigned long int loopcount = 0;
   char alphasign, alphacum;
   
   // if numloops > 0 on entry, do exactly that many iterations
   // if numloops = 0 on entry, proceed until all parities check,
   //    or until maxnumloop iterations
   if(numloops) { // set to do that many iterations
	  maxnumloop = numloops;
   }

   do {
	  loopcount++;
   	  // initialize the indexing for dealing with sparse matrices
	  bzero(na,N*sizeof(unsigned int));   // zero out this array

	  // Check node to variable node computation
	  for(m = 0; m < M; m++) {  // for each check
		 // work over nonzero elements of this row
		 l = 0;
		 n = Nm[m][l];
		 alphasign = alphantom[na[n]][n];
		 alpharow[l] = alphasign;
		 alphacum = alphasign;
		 betarow[l] = phi(betantom[na[n]][n]);
		 for(l = 1; l < Nmlen[m]; l++) {
			n = Nm[m][l];
			alphasign = alphantom[na[n]][n];
			alphacum ^= alphasign;
			alpharow[l] = alphasign;
			phibetarow[l] = phiq(betantom[na[n]][n]);
		 }
		 // compute leave-one-out min data
		 forbacksum(0,phibetarow,Nmlen[m],phioutphisum,sumall);
		 for(l = 0; l < Nmlen[m]; l++) { // for each n on this row
			n = Nm[m][l];
			alphasign = alphacum ^ alphantom[na[n]][n]; 
			phiLmton[n][na[n]++] = (1-2*alphasign)*phiq(phioutphisum[l]);
		 }
	  }  // for m
	  
      if(printstuff) printsparsetranspose("phiLmton",phiLmton);

	  // Variable node to check node computation
	  for(n = 0; n < N; n++) {  // work over each column
		 forbacksum(phiLc[n],phiLmton[n],Mnlen[n],poutsum,sumall);
		 phiLcout[n] = sumall;
		 if(phiLcout[n] >= 0) c[n] = 0; else c[n] = 1;
		 for(l = 0; l < Mnlen[n]; l++) {
			betantom[l][n] = abs(poutsum[l]);
			alphantom[l][n] = poutsum[l] < 0;
		 }
	  } // for n
	  if(printstuff) printsparse("betantom",betantom);
	  if(printstuff) printsparsechar("alphantom",alphantom);
	  if(printstuff) VECDUMP(phiLcout,N);
	  if(printstuff) VECDUMP2(c,N,int);

	  if(!numloops) { // numloops 0 on entry -- check parity
		 parity = checkparity(c);
		 if(parity) { // if parities check
			numloops = loopcount;
			return parity;
		 }
	  }
   } while(loopcount < maxnumloop);
   if(numloops) { // parity has not been previously checked
	  parity = checkparity(c);
   }
   numloops = loopcount;
   return parity;
}



void
LDPCDECODER::buildqphitable(void)
{
   Mquantphi = 2048;  // number of quantization levels
   xquantphimax = 10;
   phiquantdelta = xquantphimax/double(Mquantphi);

   qphitable = new double[Mquantphi];
   for(int i = 0; i < Mquantphi; i++) {
	  double x = i*phiquantdelta + phiquantdelta/2;
	  double phix = phi(x);
	  qphitable[i] = phix;
   }
}


double LDPCDECODER::phi(double x)
// compute the exact phi function
{
   return -log(tanh(x/2));
}

double LDPCDECODER::phiq(double x)
// compute the phi function using the interpolation table
{
   if(x > xquantphimax)
	  return qphitable[Mquantphi-1]; // return last element in table if too big
   int delta = floor(x/phiquantdelta);
   return qphitable[delta];
}

void sigtest(int sig)
{
   cout << "in sigtest sig=" << sig;
   sigflag = 1;
}

double LDPCDECODER::boxplus(double x, double y)
{
   if(x == DBL_MAX) return y;
   if(y == DBL_MAX) return x;
   return 2*atanh(tanh(x/2)*tanh(y/2));
}


double LDPCDECODER::rcbpboxplus(double a, double b)
{
   if(a == DBL_MAX) return b;
   if(b == DBL_MAX) return a;
   double N = 32;
   double delta = 0.5;
   double delta2 = 0.25;
   int r, c;
   c = std::min(floor(b/delta),(N-1));
   r = std::min(floor(a/delta),(N-1));
   if(c > r) {  // swap to make r > c
	  int rsave = r;
	  r = c;
	  c = rsave;
   }
   int d = r-c;
   int That;
   if(c > 1) {
	  That = (c-1) + (d > 2);
   }
   else {
	  That = (d > 1) & (c & 1);
   }
   double quant = That*delta + delta2;
   // printf("a=%g b=%g  r=%d  c=%d  That=%d  quant=%g  boxplus=%g\n",a,b,r,c,
   //		  That,quant,boxplus(a,b));
  return quant;
}

//---------------------------------------------------------------
#ifdef DOLPSTUFF

void vdump(vector <unsigned int> S)
{
   for(auto p = S.begin(); p != S.end(); ++p) {
	  cout << *p << " ";
   }
}

void LDPCDECODER::LPinit(void)
// initialize data structure for Linear Programming decoding
{

   int i, i1, l, m, n;
   int maxalloc = 0;
   unsigned int *LPcombos = 0;
   unsigned int *LPncomb = new unsigned int[M];
   unsigned int ncombtotal = 0;
   unsigned int totalrowwt = 0;
   char name[200];
   vector < vector < vector < unsigned int> > > E;
   vector < vector < unsigned int> > E1;

   unsigned int *ncomb;
   // E[m][combination number][element]

   unsigned int constraintmatweight = 0;
   unsigned int wtdncomb;

   for(int m = 0; m < M; m++) {
	  totalrowwt += Nmlen[m];
	  LPncomb[m] = n_even_comb(Nmlen[m],wtdncomb);
	  constraintmatweight += Nmlen[m] + wtdncomb + LPncomb[m];
// cout << "m=" << m << "  ncombs=" << LPncomb[m] << "   wtdncomb=" << wtdncomb << "  Nmlen=" <<  Nmlen[m] << "   constraintmatweight=" << constraintmatweight << endl;
	  
	  ncombtotal += LPncomb[m];
	  if(LPncomb[m] > maxalloc) {
		 maxalloc = LPncomb[m];
		 if(LPcombos) delete[] LPcombos;
		 LPcombos = new unsigned int[maxalloc + 2];
	  }
	  E1.clear();
	  for(int t = 0; t <= Nmlen[m]; t+= 2) {
		 gencomb(Nmlen[m], t, LPcombos, Nm[m], E1);
	  }
// cout << "m=" << m << "lenght(E1)=" << E1.size() << endl;
//  for(i = 0; i < E1.size(); i++) {
// 	cout << "i=" << i << "  ";
// 	vdump(E1[i]); cout << "  ";
//  }
//  cout << endl;
 
	  E.push_back(E1);
   }
// // print out the set of permutations
//  for(m = 0; m < M; m++) {
// 	for(l = 0; l < E[m].size(); l++) {
// 	   cout << "m=" << m << "  l=" << l << "   ";  
// 	   vdump(E[m][l]);
// 	   cout << endl;
// 	}
//  }
// cout << "constraintmatweight=" << constraintmatweight << endl; 
   // allocate list data for sparse constraint matrix description
   int *ia = new int[constraintmatweight+1];
   int *ja = new int [constraintmatweight+1];
   double *ar = new double [constraintmatweight+1];

   lp = glp_create_prob();
   glp_set_prob_name(lp,"LP decoder");
   glp_set_obj_dir(lp,GLP_MIN);
   glp_add_rows(lp,totalrowwt+M);
// cout << "LP rows=" << totalrowwt + M << endl;
   
   glp_add_cols(lp, N + ncombtotal);
// cout << "LP cols=" << N + ncombtotal << endl;

   // Set up constraint names and bounds
   int rowctr = 1;
   int sec = 1;  // sec = sparse element count
   //        relationships between f and w
   int start_w_offset = N;
   for(m = 0; m < M; m++) {
	  for(l = 0; l < Nmlen[m]; l++) {
		 sprintf(name,"FW:m%dl%d",m,l);
		 glp_set_row_name(lp,rowctr,name);
		 glp_set_row_bnds(lp,rowctr,GLP_FX,0,0);
		 ia[sec] = rowctr;    ja[sec] = Nm[m][l] +1;  ar[sec] = -1;
// cout << "sec=" << sec << "  ia=" << ia[sec] << "  ja=" << ja[sec] <<
// "  ar=" << ar[sec] << endl;

		 sec++;
		 int wo = 0;
		 for(i = 0; i < LPncomb[m]; i++) {
// cout << "m=" << m << "  i=" << i << "S="; vdump(E[m][i]); cout << endl;
			for(i1 = 0; i1 < E[m][i].size(); i1++) {
			   if(E[m][i][i1] == Nm[m][l]) {
				  ia[sec] = rowctr;  ja[sec] = start_w_offset + wo + 1;
				  ar[sec] = 1;
// cout << "sec=" << sec << "  ia=" << ia[sec] << "  ja=" << ja[sec] <<
// "  ar=" << ar[sec] << endl;
				  sec++;
			   }
			} // i1
			wo++;
		 } // i
		 rowctr++;
	  } // l
	  start_w_offset += LPncomb[m];
   } // m
   //         sum constraints on w
   start_w_offset = N;
// cout << "**********weight data" << endl;
   for(m = 0; m < M; m++) {
	  sprintf(name,"W:m%d",m);
	  glp_set_row_name(lp,rowctr,name);
	  glp_set_row_bnds(lp,rowctr,GLP_FX,1,1);
	  for(l = 0; l < LPncomb[m]; l++) {
		 ia[sec] = rowctr; ja[sec] = start_w_offset + l + 1;
		 ar[sec] = 1;
// cout << "sec=" << sec << "  ia=" << ia[sec] << "  ja=" << ja[sec] <<
// "  ar=" << ar[sec] << endl;
         sec++;
	  }
	  rowctr++;
	  start_w_offset += LPncomb[m];
   }
   sec--;
// cout << "sec=" << sec << endl;
   glp_load_matrix(lp,sec,ia,ja,ar);

   // Set up variable names and bounds
   int colctr = 1;
   //        code bit variables
   for(n = 0; n < N; n++) {
	  sprintf(name,"n%d",n);
	  glp_set_col_bnds(lp,colctr,GLP_DB,0,1);
	  glp_set_col_name(lp,colctr++,name);
   }
   //        even parity constraint variables
   for(m = 0; m < M; m++) {
	  for(l = 0; l < LPncomb[m]; l++) {
		 sprintf(name,"m%dS%d",m,l);
		 glp_set_col_name(lp,colctr,name);
		 glp_set_col_bnds(lp,colctr,GLP_DB,0,1);
		 colctr++;
	  }
   }

   // double *lambda = new double[N];
   // double *f = new double[N];
   // lambda[0] = 1;  lambda[1] = 1; lambda[2] = 1; lambda[3] = -1;
   // lambda[4] = -1;  lambda[5] = -1; lambda[6] = -10;

   // for(n = 0; n < N; n++) {
   // 	  glp_set_obj_coef(lp,n +1, lambda[n]);
   // }
   for(n = N; n < N + ncombtotal; n++) {
	  glp_set_obj_coef(lp, n +1, 0);
   }
   // It appears that once these coefficients are set, they
   // don't need to be set again for another problem.
   // So only the first N coefficients need to be set for
   // a decoding problem

   // glp_simplex(lp,NULL);
   // double z;
   // z = glp_get_obj_val(lp);
   // for(n = 0; n < N; n++) {
   // 	  f[n] = glp_get_col_prim(lp,n +1);
   // }
   // cout << "z=" << z << "  ";
   // VECDUMP(f,N);

   // // explore the interface a bit
   // cout << "objective coefficients: ";
   // for(n = 0; n < N + ncombtotal; n++) {
   // 	  cout << glp_get_obj_coef(lp,n +1) << "  ";
   // }
   // cout << endl;
   
   // clean up
   // glp_delete_prob(lp);

   
   // set parameters so that messaging is off
   glp_init_smcp(&glp_param);
   glp_param.msg_lev = GLP_MSG_OFF;

   delete [] LPcombos;
   delete [] LPncomb;
   delete[] ia;
   delete [] ja;
   delete [] ar;
   // delete [] lambda;
   // delete [] f;
}

unsigned int LDPCDECODER::n_even_comb(unsigned int n,
									  unsigned int& weightedncomb)
{
   unsigned int ncomb = 0;		// maximum number of combinations
   unsigned int ncomb1;
   weightedncomb = 0;
   for(unsigned int t = 0; t <= n; t+= 2) {
	  ncomb1 = nchoosek(maxrowwt,t);
	  weightedncomb += ncomb1 * t;
	  ncomb += ncomb1; 
   }
   return ncomb;
}

unsigned int LDPCDECODER::nchoosek(unsigned int n, unsigned int k)
// compute n choose k = n!/(k!(n-k)!))
// (won't work with numbers that are too large
{
   unsigned int nfact;
   unsigned int kfact;
   unsigned int nkfact;
   nfact = factorial(n);
   kfact = factorial(k);
   nkfact = factorial(n-k);
   return nfact/(kfact*nkfact);
}

unsigned int LDPCDECODER::factorial(unsigned int n)
{
   unsigned int p = 1;
   for(unsigned int i = 2; i <= n; i++) {
	  p *= i;
   }
   return p;
}


void LDPCDECODER::gencomb(int n, int t, unsigned int *c, int *Nm,
						   vector <vector <unsigned int>>& E1)
// generate all t-combinations of the n numbers {0,...,n-1}
// each combination if c_t, ... c_2, c_1
// Algorithm 7.2.1.3L of Knuth (V4 Fascicle 3)
// additional variables c_{t+1} and c_{t+2} are used as sentinels
{
   int j;
   vector< unsigned int> S;
   // L1: initialize
   for(j = 0; j < t; j++) {
	  c[j] = j;
   }
   c[t] = n;
   c[t+1] = 0;
   while(1) {
	  // L2: visit the combination c_t ... c_2 c_1
	  S.clear();
	  for(int i = 0; i < t; i++) {
		 // cout << c[i] << " ";
		 S.push_back(Nm[c[i]]);
	  }
// cout << "S=";
// for(auto p = S.begin(); p != S.end(); ++p) {
// 	 cout << *p << " ";
// }
// cout << endl;
// 	  cout << endl;
	  E1.push_back(S);
	  // L3: Find j
	  j = 0;
	  while(c[j]+1 == c[j+1]) {
		 c[j] = j;
		 j++;
	  }
	  if(j >= t) break;
	  c[j] = c[j] + 1;
   }
   // cout << "E1=" << E1 << endl;
   

}

// ------------------------------------------------------------
// Stuff for LP decoding
void LDPCDECODER::LPdecodeinit(double *y)
{
   int n,m;
   double Lc;
   Lc = 2*a/sigma2;
   //  set up data to test results by hand for testing
   // Lc = 1;
   // y[0] = 1; y[1] = 1; y[2] = 1;
   // y[3] = -1; y[4] = -1; y[5] = -1;
   // y[6] = 1;
   
   for(n = 0; n < N; n++) {
	  glp_set_obj_coef(lp,n +1, Lc*y[n]);
   }
}

int LDPCDECODER::LPdecode(unsigned char *c, unsigned long int maxnumloop,
							unsigned long int &numloops)
{
   int n;
   int parity;

   glp_simplex(lp,&glp_param);
   double z;
   double f;
   z = glp_get_obj_val(lp);
// cout << "LP: z = " << z << "  x=";
   for(n = 0; n < N; n++) {
	  f = glp_get_col_prim(lp,n +1);
// cout << f << "  ";
	  if(f > 0.5) c[n] = 1;
	  else c[n] = 0;
   }
// cout << endl;
   parity = checkparity(c);
   numloops = 1;   // one LP simplex algorithm
   return parity;
}



#endif  // DOLPSTUFF


/*
Local Variables:
compile-command: "g++ -c -g -std=c++11 ldpcdecoder.cc -Wno-deprecated"
End:
*/

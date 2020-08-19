//
//  Program: polarcode -- functions for polar codes
//       encoding, SC decoding
//       systematic encoding, systematic decoding
//
//
//  Todd K. Moon
//  Utah State University
//
//  Date: Started 30 June 2018//
//


// A note on structure and design:  This file incorporates several encoders and
// decoders.
// Encoding main functions:
//           encode(  )              (nonsystematic encoding)
//           encodesyst(  )          (systematic encoding)
//           encodewithCRC(  )       (nonsystematic encoding with CRC)
//           encodesystwithCRC(   )  (systematic encoding with CRC)
//
// Decoding main functions:
//           SCdecode(   )           (successive cancellation decoder.
//                                   (does not return decoded codeword)
//           decodesyst(  )          (decode systematically encoded codewords)
//           listdecode(   )         (list decoder using likelihoods)
//           listdecodeLLR(  )       (list decoder using log likelihoods)
//  Each of these decoding functions is able to do CRC decoding,
//     when properly called.
//  In terms of design, it might have been better to do this
//  using a class hierarchy than a single class, but there is some
//  efficient memory use in the way this was set up.  You just have to
//  be careful with how you call things.

// There is a single constructor for the class.  There are arguments
//  that indicate to it how to allocate memory.
//
// polarcode(int N_in, int Kp_in,
//           encodertype enctype_in,
//   		 decodertype dectype_in
//			 int L_in = 1, int crcsize_in = 0);
//
// enctype_in: tells what encoder type may be used, and allocates memory
//    accordingly.  Options are: nonsystematicenc, systematicenc, withcrcenc.
//  These can be OR-ed together if multiple capabilities are desired with a single
//  object.  For example, you can use
//      enctype_in = nonsystematicnec | systematicenc | withcrcenc
//  to create an object with memory allocated for both systematic and
//  nonsystematic encoding, with CRC encoding as well.
//
//  Note that enctype_in is used in the construction process, and you determine
//    which encoder is used by the function you call.
//
// dectype_in: tells what decoder may be used, and allocates memory
//   accordingly.  Options are: SCdec, listdec, listdec, withcrcdec.
//  These can be OR-ed together if multiple capabilities are desired with a single
//  object.  For example,
//       dectype_in = SCdec | listdec | listllrdec | withcrcenc
//
// Note that detype_in is used in the consctruction process, and you determine
//   which decoder is used by the function you call.
// The decoders also accept an argument thisenctype indicating what kind of decoding is done.
// This can be used to tell the list-based decoders if CRC enoding is used.
//  (this argument has no effect on the SCdecode function.)
// The decoders return the decoded message (Kp bits long) and the
// decoded codeword (N bits long) via arguments passsed by reference.
//  The SCdecode function only returns the decoded message, not the codeword,
//  due to the structure of the decoding algorithm.  The codeword can be obtained in
//  this case by re-encoding the message bits.
//  The decoding functions return a bool indicating whether the codeword is returned.
//  Thus, SCdecode returns false, while the list-based decoders return true.
//
// Designing polar codes:
//   designBhatt( ) -- design using iterated Battacharrya transformations
//      z --> z^2, 2z - z^2
//    
//   designBhatt2( ) -- design using iterated Battacharrya transformations
//      z --> z^2, 2z - z^2
//    works in the bit-reversed space, so no reversal is necessary
//
//   designMC( ) -- design by simulating bits separately, 
//      then choosing the most reliable bits




#include "polarcode.h"
#include <iostream>
#include <cstring>  // memset
#include <limits>   // max()
using namespace std;

#include "printstuff.cc"


polarcode::polarcode(int N_in, int Kp_in,
					 encodertype enctype_in,
					 decodertype dectype_in, int L_in, int crcsize_in)
// Constructor:
// Kp: size of code without the CRC (if any).  K = Kp + crcsize
// If L_in > 1, then list decoding is assumed.
// If crcsize_in != 0, then crc decoding is assumed (relevant only for list decoding)
{

   // check some consistency
   if(!(enctype_in & (nonsystematicenc | systematicenc))) {
	  cerr << "Must select either nonsystematic or systematic encoding" << endl;
	  exit(-1);
   }
   if((!(dectype_in & (listdec | listllrdec)) && L_in > 1)) {
   	  cerr << "Error: SC decoder with L=" << L_in << endl;
   	  exit(-1);
   }
   if(crcsize_in && crcsize_in != 16) {
	  cerr << "only 16-bit CRC presently implemented" << endl;
	  exit(-1);
   }
   if(crcsize_in && !(enctype_in & withcrcenc)) {
	  cerr << "CRC requested without crc decoding" << endl;
	  exit(-1);
   }
   if(L_in == 0 && (dectype_in & (listdec | listllrdec))) {
	  cerr << "List decoding with list of size 0" << endl;
	  exit(-1);
   }
   if(L_in && !(dectype_in & (listdec | listllrdec))) {
	  cerr << "List size=" << L_in << "  but no list decoding requested" << endl;
	  exit(-1);
   }
   
   N = N_in;
   K = Kp_in + crcsize_in;
   n = int(log2(N) + 0.1);		// number of levels in tree
   L = L_in;
   crcsize = crcsize_in;
   do_message = false;
   if(crcsize) do_message = true;  // relevant only for list decoding

   enctype = enctype_in;
   dectype = dectype_in;

   // Data relevant to all encoders/decoders
   polarcodedesign = new SIGNEDBITTYPE[N];
   bitreverser = new int[N]; // array of bit-reversed data

   bitsout = new BITTYPE[N];
   messagebitsout = new BITTYPE[N];
   uwork = new BITTYPE[N];		// working space for encoder
   uworksave = new BITTYPE[N];	// save input message for test

   // set indicators of what is not allocated
   // SC decoder
   level1array = 0;
   level0array = 0;
   LLR = 0;
   LLRfinal = 0;
   BITS = 0;
   Bsys = 0;
   // SCL decoder
   pathIndexToArrayIndex = 0;
   arrayReferenceCount = 0;
   pathmetric = 0;
   inactiveArrayIndices = 0;
   activePath = 0;
   arrayPointer_P = 0;
   arrayPointer_C = 0;
   arrayPointer_llr = 0;
   probForks = 0;
   contForks = 0;
   idx = 0;
   arrayPointer_message = 0;
   
   int level1, level0;
   for(int i = 0; i < N; i++) {
	  polarcodedesign[i] = -1; 	// default to not frozen bits everywhere
	  bitreverser[i] = bitreverse(i,n,level1,level0);
   }

   if(dectype & SCdec) {			// successive cancellation decoding
	  level1array = new int[N];	// (actually, used for SC decoder)
	  level0array = new int[N];
	  for(int i = 0; i < N; i++) {
		 bitreverse(i,n,level1,level0);
		 level1array[i] = level1;
		 level0array[i] = level0;
	  }
	  LLR = new double[2*N-1];
	  LLRfinal = new double[N];
	  CALLOCMATRIX(BITS,BITTYPE,2,N-1);
   }

   setdoprintvalue(1);

   if(enctype & systematicenc){ 
	  CALLOCMATRIX(Bsys,BITTYPE,n+1,N);
	  // visually, this should be N x (n+1)  (the shape of the
	  //   encoding graph.
	  //  However, I transpose this to provide some arrays
	  //    of length N to be used as workspace.
   }

   if(dectype & (listdec | listllrdec)) {
	  // if a list decoder
	  buildListDecodeDataStructures();
   }

   // set the CRC code information (no extra space)
   crcp = 16;
   crcg = 0x18005;//CRC divisor: 1 1000 0000 0000 0101 = x^16+x^15+x^2+1
   crcmask = (1<<crcp) - 1;

}

polarcode::~polarcode()
{
   delete[] polarcodedesign;
   delete [] bitreverser;
   delete [] bitsout;
   delete [] messagebitsout;
   delete [] uwork;
   delete [] uworksave;
   if(probForks) { FREEMATRIX(probForks); }
   if(contForks){ FREEMATRIX(contForks); }
   if(idx) delete[] idx;


   if(level1array) delete[] level1array;
   if(level0array) delete[] level0array;
   if(LLR) delete[] LLR;
   if(LLRfinal) delete[] LLRfinal;
   
   if(BITS) { FREEMATRIX(BITS);}
   if(Bsys) { FREEMATRIX(Bsys); }
   if(pathIndexToArrayIndex) { FREEMATRIX(pathIndexToArrayIndex); }
   if(arrayReferenceCount) { FREEMATRIX(arrayReferenceCount); }
   if(inactiveArrayIndices) delete[] inactiveArrayIndices;
   if(activePath) delete[] activePath;
   if(pathmetric) delete[] pathmetric;
   if(arrayPointer_P) {
	  for(int lambda = 0; lambda <= n; lambda++) {
		 for(int s = 0; s < L; s++) {
			FREEMATRIX(arrayPointer_P[lambda][s]);
		 }
	  }
	  FREEMATRIX(arrayPointer_P);
   }
   if(arrayPointer_C) {
	  for(int lambda = 0; lambda <= n; lambda++) {
		 for(int s = 0; s < L; s++) {
			FREEMATRIX(arrayPointer_C[lambda][s]);
		 }
	  }
	  FREEMATRIX(arrayPointer_C);
   }
   if(arrayPointer_llr) {
	  for(int lambda = 0; lambda <= n; lambda++) {
		 for(int s = 0; s < L; s++) {
			delete[] arrayPointer_llr[lambda][s];
		 }
	  }
	  FREEMATRIX(arrayPointer_llr);
   }
   if(arrayPointer_message) { FREEMATRIX(arrayPointer_message); }
}


bool polarcode::SCdecode(double *y, BITTYPE* & decout,
							 BITTYPE* & u2,encodertype no_op)
// Successive cancellation decoder, following Tal/Vardy Trans IT, v. 61
// no. 5, May 2015, pp 2213--2225.
//
// returns false, indicating that the codeword is NOT returned (decout untouched)
// the K-bit message is in messagebitsout
// encodertype no_op -- a do-nothing argument, use to keep the
//    decoder argument lists the same
{
   int level1, level0;
   int i,j;			// bit index, and bit-reversed version
   
   computeLLR(y, LLR + N-1);	// write starting at location LLR[N-1]
   for(j = 0; j < N; j++) {   // for each bit in order
	  i = bitreverser[j];
	  level1 = level1array[j];
	  level0 = level0array[j];
	  updateLLR(i,level1);	// process the LLRs up the tree
	  LLRfinal[i] = LLR[0];	// save the LLR
	  LLRtobit(j,i);				// convert LLR to bit
	  updateBITS(bitsout[j],i,level0); // push the bits back down the tree
   }
   for(i = 0, j = 0; i < N; i++) {
	  if(polarcodedesign[i] == -1) {
		 messagebitsout[j++] = bitsout[i];
	  }
   }
   u2 = messagebitsout;
   return false;
}

void polarcode::updateLLR(int i, int level1)
// part of successive cancellation decoding
{
   int start, end, nextlevel, lastlevel;

   if(i == 0) {
	  nextlevel = n;
   }
   else {
	  lastlevel = level1;
	  start = (1<<(lastlevel-1)) - 1;
	  end = (1 << lastlevel) - 1;
	  // At the deepest level of the tree, do the LLRlower 't' computations
	  for(int index = start; index < end; index++) {
		 LLR[index] = LLRlower(BITS[0][index], LLR[end + 2*(index - start)],
							   LLR[end + 2*(index - start) + 1]);
	  }
	  nextlevel = level1 - 1;
   }
   // At the other levels of the tree, do the LLRupper 's' computations
   for(int level = nextlevel; level >= 1; level--) {
	  start = (1<<(level-1)) - 1;
	  end = (1 << level) - 1;
	  for(int index = start; index < end; index++) {
		 LLR[index] = LLRupper(LLR[end + 2*(index - start)],
							   LLR[end + 2*(index - start) + 1]);
	  }
   }
}

void polarcode::updateBITS(BITTYPE latestbit, int i, int level0)
// part of successive cancellation decoding
{
   int start, end, level, index, lastlevel;
   
   if(i == N-1) {
	  return;		// no processing to do on the last
   }
   else if(i < N/2) {
	  BITS[0][0] = latestbit;
   }
   else {
	  // Compute and propagate bits in the "lower" part
	  lastlevel = level0;
	  BITS[1][0] = latestbit;
	  for(level = 1; level <= lastlevel - 2; level++) {
		 start = (1<<(level-1)) - 1;
		 end = (1 << level) - 1;
		 for(int index = start; index < end; index++) {
			BITS[1][end + 2*(index - start)] = (BITS[0][index] + BITS[1][index]) % 2;
			BITS[1][end + 2*(index - start) + 1] = BITS[1][index];
		 }
	  }
	  // Compute and propagate bits in the "upper" part
	  level = lastlevel - 1;
	  start = (1<<(level-1)) - 1;
	  end = (1 << level) - 1;
	  for(index = start; index < end; index++) {
		 BITS[0][end + 2*(index-start)] = (BITS[0][index]+BITS[1][index]) % 2;
		 BITS[0][end + 2*(index-start) +1] = BITS[1][index];
	  }
   }
}

	  
double polarcode::LLRupper(double LLR1, double LLR2)
// part of successive cancellation decoding
{
   // Use of the logdomain_sum here makes it so that the argument
   // of every exponential exp(x) is negative.  This stabilizes
   // the computations, and also suggests ways of approximating the
   // terms.

   return logdomain_sum(LLR1 + LLR2, 0) - logdomain_sum(LLR1, LLR2);
}


double polarcode::LLRlower(BITTYPE bit, double LLR1, double LLR2)
// part of successive cancellation decoding
{
   if(bit == 0) {
	  return LLR1 + LLR2;
   }
   else {
	  return LLR2 - LLR1;
   }
}

double polarcode::logdomain_sum(double x, double y)
{
   if(y > x) {
	  return y + log(1 + exp(x - y));
   }
   else {
	  return x + log(1 + exp(y - x));
   }
}

double polarcode::logdomain_diff(double x, double y)
// x must be greater than y
{
   return x + log(1 - exp(y - x));
}


void polarcode::LLRtobit(int j, int i)
// convert Log likelihood ratios to bits, taking the polar code design into account
{
   if(polarcodedesign[j] == -1) {
	  if(LLR[0]  > 0) {
		 bitsout[j] = 0;
	  }
	  else {
		 bitsout[j] = 1;
	  }
   }
   else {
   	  bitsout[j] = polarcodedesign[j];
   }
}

  
   
int polarcode::bitreverse(int x, int n, int& level1,
						   int& level0)
// return the bit-reversed version of x.  Also, level1, which is
// the index of the first 1, and level0, which is the index of the first 0
{
   // level1:
   // level0: the index of the first zero
   int i;
   unsigned char bits[64];
   int mask = 1;
   int out = 0;
   bool onefound = false;
   bool zerofound = false;
   level1 = n;
   level0 = 0;
   
   for(int i = 0; i < n; i++) {
	  if(mask & x) {
		 bits[i] = 1;
	  }
	  else {
		 bits[i] = 0;
	  }
	  mask <<= 1;
   }
   mask >>= 1;
   for(i = 0; i < n; i++) {  // i == 0 is MSB of the output
	  if(bits[i]) {
		 out += mask;
		 // also find the first index of 1
		 if(!onefound) {
			level1 = i+1;
			onefound = true;
		 }
	  }
	  else {
		 if(!zerofound) {
			level0 = i + 1;
			zerofound = true;
		 }
	  }
	  
	  mask >>= 1;
   }
   return out;
}


BITTYPE* polarcode::encode(BITTYPE* u)
// Encode data nonsystematically.  Returns a pointer to the codeword
// Does not do CRC
{

   // Merge u in with the frozen data and bit-reverse permute
   for(int i = 0, j=0; i < N; i++) {
	  if(polarcodedesign[i] == -1) { // not a frozen bit
		 uwork[bitreverser[i]] = u[j++];
	  }
	  else {					// a frozen bit
		 uwork[bitreverser[i]]= polarcodedesign[i];
	  }
   }
   for(int i = 0; i < N; i++) { // save the input message (for debug/test)
	  uworksave[i] = uwork[bitreverser[i]];
   }
   vectFn(uwork);  
   return uwork;
}

void polarcode::vectFn(BITTYPE *uwork)
// multiply uwork by F kron^n (n-fold kronecker product of F)
// product is computed inplace (so uwork is modified)
{
   int skipr = 1;  // offset on RHS of adder equation.  Also number of addergroups
   int stepl = 2;  // left index steps
   int ninaddergroup = N/2;  // number of adders in group in current level
   int lindexstart = 0;  // starting index of addergroup

   // in-place computation (answer is in uwork)
   for(int n1 = 0; n1 < n; n1++) { // for each level
	  for(int i2 = 0; i2 < skipr; i2++) {
		 int lindx = lindexstart;
		 for(int i1 = 0; i1 < ninaddergroup; i1++) {
			uwork[lindx] = uwork[lindx] ^ uwork[lindx+skipr];  // mod 2 sum
			lindx += stepl;
		 } // end for i1
		 lindexstart++;
	  }  // end for i2
	  skipr <<= 1;
	  stepl <<= 1;
	  ninaddergroup >>= 1;
	  lindexstart = 0;
   } // end for n1
}
   
BITTYPE* polarcode::unencode(BITTYPE *x)
// given an ostensive codeword x, return the corresponding message u
// Performs the F^n transformation, followed by bit reverse permutation
{

   vectFn(x);    // run through unitary part
   for(int i = 0; i < N; i++) {
	  uworksave[i] = x[bitreverser[i]];
   }
   
   for(int i = 0, j=0; i < N; i++) {
	  if(polarcodedesign[i] == -1)
		 messagebitsout[j++] = uworksave[i];
   }
   return messagebitsout;
}

void polarcode::designBhatt(double designSNRdB)
// Design using the Battacharyya parameters Z --> 2Z - Z^2    Z --> Z^2
// this version works in the bit-permuted index space,
// so it is reversed before sorting.  (Following Vangala, Viterbo & Hong)
{
   double *z = new double[N];	// the Bhattacharyya parameters computed
   
   double tz;					// temporary z
   double l2 = log(2.0);
   int *idx = new int[N];
   void sort2double_int(int num_elements, double *array, int *array2);

   double *zreverse = new double[N];

   for(int i = 0; i < N; i++) {
	  idx[i] = i;
	  polarcodedesign[i] = 0;
   }
   
   double designSNR = pow(10,designSNRdB/10);
   z[0] = -designSNR;

   int B = 2;   // number of nodes at first level
   while(B <= N) {
	  int B2 = B/2;
	  for(int j = 0; j < B2; j++) {
		 tz = z[j];
		 z[j] = logdomain_diff( l2 + tz, 2*tz);
		 z[B2 + j] = 2*tz;
	  }
	  B <<= 1;
   }

   // bit reverse assignments
   for(int i = 0; i < N; i++) {
	  zreverse[i] = z[bitreverser[i]];
   }
   for(int i = 0; i < N; i++) {
	  z[i] = zreverse[i];
   }

   sort2double_int(N,z,idx);
   
   for(int i = 0; i < K; i++) {
	  polarcodedesign[idx[i]] = -1;
   }

   delete [] zreverse;
   delete[] z;
   delete[] idx;
}

void polarcode::designBhatt2(double designSNRdB)
// Design using the Battacharyya parameters Z --> 2Z - Z^2    Z --> Z^2
// This version works in the original bit order domain.
{
   double *z = new double[N];	// the Bhattacharyya parameters computed
   
   double tz;					// temporary z
   double l2 = log(2.0);
   int *idx = new int[N];
   void sort2double_int(int num_elements, double *array, int *array2);

   for(int i = 0; i < N; i++) {
	  idx[i] = i;
	  polarcodedesign[i] = 0;
   }
   
   double designSNR = pow(10,designSNRdB/10);
   z[0] = -designSNR;

   int nstep = 1;				// number of steps at each level
   int lskip = N/2;				// skip on index on LHS
   int rskip = N;				// skip on index on RHS
   int rindex;

   for(int j = 0; j < n; j++) {		// for each level
	  rindex = 0;
	  for(int t = 0; t < nstep; t++) {
		 tz = z[rindex];
		 z[rindex] = logdomain_diff( l2 + tz, 2*tz);
		 z[rindex + lskip] = 2*tz;
		 rindex += rskip;
	  }
	  lskip >>= 1;
	  rskip >>= 1;
	  nstep <<= 1;
   }


   sort2double_int(N,z,idx);
   printdoublevec("zBhatt",z,N);

   
   for(int i = 0; i < K; i++) {
	  polarcodedesign[idx[i]] = -1;
   }

   delete[] z;
   delete[] idx;
}

void polarcode::designMC(double designSNRdB, int M)
// Design using Monte Carlo method.  M is the number of bits to estimate bit error
{
   double *y = new double[N];	// simulated channel received value
   int *idx = new int[N];
   int *MCbiterrcount = new int[N];
   void sort2int_int(int num_elements, int *array, int *array2);
   double gran(void);

   for(int i = 0; i < N; i++) {
	  idx[i] = i;
	  polarcodedesign[i] = 0;  // freeze all the bits
	  MCbiterrcount[i] = 0;
   }
   double R = double(K)/double(N);
   double Ec = 1;				// coded BPSK energy
   double Ecsqrt = sqrt(Ec);
   double Eb = Ec/R;			// energy per bit
   double EbN0 = pow(10.0, designSNRdB/10);
   double N0, sigma2, sigma;
   BITTYPE *u_decoded;
   BITTYPE *cw;

   N0 = Eb/(EbN0); sigma2 = N0/2;  sigma = sqrt(sigma2);
   setchannelfact(-4*sqrt(Ec)/N0);

   for(int bitno = 0; bitno < N; bitno++) { // simulate each bit separately,
	                                        // with the other bits frozen
	  if(bitno) polarcodedesign[bitno-1] = 0; // freeze previous bit
	  polarcodedesign[bitno] = -1;			  // unfreeze current bit
	  for(int m = 0; m < M; m++) {			  // simulate M bits
		 for(int i = 0; i < N; i++) y[i] = -Ecsqrt + sigma*gran(); // all-zero 
		 SCdecode(y,cw,u_decoded,0);
		 if(u_decoded[m] == 1) { // decoded incorrectly
			MCbiterrcount[bitno]++;
		 }
	  }
   }
   
   sort2int_int(N,MCbiterrcount,idx);

   printintvec("MCbiterrcount",MCbiterrcount,N);
   
   for(int i = 0; i < K; i++) {
	  polarcodedesign[idx[i]] = -1;
   }

   delete[] MCbiterrcount;
   delete[] y;
   delete[] idx;
}

void polarcode::setdesign(SIGNEDBITTYPE *polarcodedesign_in)
// Set the design of this object according to the design passed in
{
   for(int i = 0; i < N; i++) {
	  polarcodedesign[i] = polarcodedesign_in[i];
   }
}


BITTYPE *polarcode::encodesyst(BITTYPE *u)
// Encoding systematically.  Returns a pointer to the codeword
// Sets up the system of equations, then solves by passing bits across
// the encoding strucure.  Does the bit reverse permutation.
{
   // fill column 0 of Bsys with 0 for indices not in A
   // fill column n of Bsys with codebits for indices in A
   memset(Bsys[0],0,sizeof(BITTYPE)*N*(n+1));
   for(int i = 0, j = 0; i < N; i++) {
	  if(polarcodedesign[i] == -1) { // not frozen
		 Bsys[n][i] = u[j++];
	  }
	  else { // frozen
		 Bsys[0][i] = 0;
	  }
   }
   int s, delta, t, l, bl, kappa;
   for(int i = N-1; i >= 0; i--) {
	  if(polarcodedesign[i] == -1) { // not frozen
		 s = n-1;
		 delta = -1;
	  }
	  else { 					// frozen
		 s = 1;
		 delta = 1;
	  }
	  t = s;
	  for(int j = 0; j < n; j++) { // for each level needing computation
		 l = MIN(t, t-delta);
		 bl = (i & (1 << (n-l-1)));
		 if(bl == 0) { 			// a line with an adder at this level
			kappa = (1 << (n-1-l));
			Bsys[t][i] = Bsys[t - delta][i] ^ Bsys[t - delta][i + kappa];
		 }
		 else { 				//  a line with no adder at this level -- copy
			Bsys[t][i] = Bsys[t - delta][i];
		 }
		 t += delta;   // t = s + j*delta
	  } // for j
	  uwork[bitreverser[i]] = Bsys[n][i];  // copy with bit reverse permutation
   } // for i

   return uwork;
}

bool polarcode::decodesyst(double *y, BITTYPE* &codewordout,
							   BITTYPE* &u2, encodertype no_op)
// Decodes a systematically-encoded codeword using successive cancellation.
// The codeword is returned in the codeword out argument.
// The information bits are returned in the u2 argument.
// no_op is a do_nothing argument, used so all decoders have the same
// argument lists.
// returns true, indicating that the codeword is computed.
{
   BITTYPE *zsyst; // first stage of decoding
   BITTYPE *cw;
   SCdecode(y,cw,zsyst,no_op);
   vectFn(zsyst);
   // extract the message bits, and bit permute to finish re-encoding
   BITTYPE *work = Bsys[0];  // grab some space to bit reverse swap and return
   for(int i =0, j = 0; i < N; i++) {
	  if(polarcodedesign[i] == -1) {
		 messagebitsout[j++] = zsyst[i];
	  }
	  work[i] = zsyst[bitreverser[i]];
   }
   u2 = messagebitsout;
   // return work;  // NOTE:  The first row of Bsys is used to hold this
   //               //  decoded codeword
   return true;
}

// ----------------------------------------------------------------------
// Stuff for list decoding

void polarcode::buildListDecodeDataStructures()
// Tal/Vardy Algorithm 5
//
// Build the data structures for LR list decoding (performed once)
{
   inactiveArrayIndices = new polarstack[n+1];
   inactivePathIndices.initspace(L);
   for(int i = 0; i <= n; i++) {
   	  inactiveArrayIndices[i].initspace(L);
   }
   // Set up data structures
   activePath = new bool[L];
   CALLOCMATRIX(arrayPointer_C,BITTYPE**, n+1,L);
   if(dectype & listdec) {
	  CALLOCMATRIX(arrayPointer_P,double **,n+1,L);
   }
   if(dectype & listllrdec) {
	  CALLOCMATRIX(arrayPointer_llr,double *,n+1,L);
	  pathmetric = new double[L];
   }
   CALLOCMATRIX(pathIndexToArrayIndex,int,n+1,L);
   CALLOCMATRIX(arrayReferenceCount,int,n+1,L);
   for(int lambda = 0; lambda <= n; lambda++) {
	  for(int s = 0; s < L; s++) {
		 if(dectype & listdec) {
			CALLOCMATRIX(arrayPointer_P[lambda][s],double,1<<(n-lambda),2);
		 }
		 if(dectype & listllrdec) {
			arrayPointer_llr[lambda][s] = new double[1<<(n-lambda)];
		 }
		 CALLOCMATRIX(arrayPointer_C[lambda][s],BITTYPE,1<<(n-lambda),2);
	  }
   }

   CALLOCMATRIX(probForks,double,L,2);  // for Alg. 13
   CALLOCMATRIX(contForks,bool,L,2);  // for Alg. 13
   idx = new int[2*L];                // for Alg. 13

   if(do_message) {
	  CALLOCMATRIX(arrayPointer_message, BITTYPE,L,N);
	  // Holds the message bits
   }
}



void polarcode::initializeListDecodeDataStructures()
//
// Initialize the data structures for LR list decoding (performed for each
// decoded codeword)
{

   for(int i = 0; i <= n; i++) {
	  inactiveArrayIndices[i].clearstack();
   }
   inactivePathIndices.clearstack();
   
   for(int lambda = 0; lambda <= n; lambda++) {
	  for(int s = L-1; s >= 0; s--) {
		 arrayReferenceCount[lambda][s] = 0;
		 inactiveArrayIndices[lambda].push(s);
	  }
   }
   for(int l = L-1; l >= 0; l--) {
   	  activePath[l] = false;
	  inactivePathIndices.push(l);
	  if(dectype & listllrdec) {
		 pathmetric[l] = 0;
	  }
   }
}


bool polarcode::listdecode(double *y,BITTYPE* &decout, BITTYPE* &info,
						   encodertype thisenctype)
// Tal/Vardy Algorithm 12
//
// List decoding using likelihoods.
//
// returns true, indicating that the decoded codeword is placed in decout
//
//  If systematic CRC encoding is used, then this decoder must be
//    informed that thisenctype = systematicenc.
{

   initializeListDecodeDataStructures();    // clear memory
   int l = assignInitialPath();

cout << "initial:" << endl;
dumpdat();
   double **P0 = getArrayPointer_P(0,l);
   for(int beta = 0; beta < N; beta++) {
	  double p0, p1, pm;
	  p0 = W(y[beta],0);
	   p1 = W(y[beta],1);
	  pm = MAX(p0,p1);
	  P0[beta][0] = p0/pm;
	  P0[beta][1] = p1/pm;
   } // for beta
   for(int phi = 0; phi < N; phi++) {
	  recursivelyCalcP(n,phi);
	  if(polarcodedesign[phi] != -1) {  // u[phi] is frozen
		 for(int lp = 0; lp < L; lp++) {
			if(activePath[lp] == false) continue;
			BITTYPE ** Cm = getArrayPointer_C(n,lp);
			Cm[0][phi & 1] = polarcodedesign[phi];
		 }
		 if(do_message) {
			arrayPointer_message[l][phi] = polarcodedesign[phi];
		 }
	  }
	  else {
		 continuePaths_UnFrozenBit(phi);
	  }
	  if(phi & 1) {
		 recursivelyUpdateC(n,phi);
	  }
cout << "phi=" << phi << endl;
dumpdat();
   }
   l = findMostProbablePath( crcsize != 0, thisenctype);
   BITTYPE ** C0 = getArrayPointer_C(0,l);
   for(int beta = 0; beta < N; beta++) {
	  decout[beta] = C0[beta][0];
   }

   if(do_message) {
	  if(thisenctype == nonsystematicenc) {
		 extractmessage(messagebitsout,arrayPointer_message[l]);
	  }
	  else if(thisenctype == systematicenc) {
		 extractmessagefromsyst(messagebitsout,decout);
	  }
	  else {
		 cerr << "thisenctype needs to be set" << endl;
	  }
	  
   }
   else {  // take the codeword and unencode it to get the message bits
	  for(int i = 0; i < N; i++) {
		 messagebitsout[i] = decout[i];   // leave decout untouched
	  }
	  if(thisenctype == nonsystematicenc) {
		 unencode(messagebitsout);
	  }
	  else if(thisenctype == systematicenc) {	
		 extractmessagefromsyst(messagebitsout,decout);
	  }
	  else {
		 cerr << "thisenctype needs to be set" << endl;
	  }
   }
   info = messagebitsout;
   return true;
}

bool polarcode::listdecodeLLR(double *y, BITTYPE* &decout, BITTYPE* &info,
							  encodertype thisenctype)
// Tal/Vardy Algorithm 12, modified for LLR decoding
//
// List decoding using log likelihood ratios.
//
// returns true, indicating that the decoded codeword is placed in decout
//
//  If systematic CRC encoding is used, then this decoder must be
//    informed that thisenctype = systematicenc.
{

   initializeListDecodeDataStructures();    // clear memory
   int l = assignInitialPath();
   double *P0 = getArrayPointer_llr(0,l);
   for(int beta = 0; beta < N; beta++) {
	  double p0, p1, pm;
	  p0 = W(y[beta],0);
	  p1 = W(y[beta],1);
	  P0[beta] = log(p0/p1);
   } // for beta
   for(int phi = 0; phi < N; phi++) {
	  recursivelyCalcllr(n,phi);
	  if(polarcodedesign[phi] != -1) {  // u[phi] is frozen
		 // continuePaths_FrozenBit
		 for(int lp = 0; lp < L; lp++) {
			if(activePath[lp] == false) continue;
			BITTYPE ** Cm = getArrayPointer_Cllr(n,lp);
			double *llr_p = getArrayPointer_llr(n,lp);
			Cm[0][phi & 1] = polarcodedesign[phi];
			pathmetric[lp] += log(1 + exp(-llr_p[0]));
		 }
		 if(do_message) {
			arrayPointer_message[l][phi] = polarcodedesign[phi];
		 }
	  }
	  else {
		 continuePaths_UnFrozenBitllr(phi);
	  }
	  if(phi & 1) {
		 recursivelyUpdateCllr(n,phi);
	  }
   } // for phi
   l = findMostProbablePathllr(crcsize != 0, thisenctype);
   BITTYPE ** C0 = getArrayPointer_Cllr(0,l);
   for(int beta = 0; beta < N; beta++) {
	  decout[beta] = C0[beta][0];
   }


   if(do_message) {
	  if(thisenctype == nonsystematicenc) {
		 extractmessage(messagebitsout,arrayPointer_message[l]);
	  }
	  else if(thisenctype == systematicenc) {
		 extractmessagefromsyst(messagebitsout,decout);
	  }
	  else {
		 cerr << "thisenctype needs to be set" << endl;
	  }
	  
   }
   else {  // take the codeword and unencode it to get the message bits
	  for(int i = 0; i < N; i++) {
		 messagebitsout[i] = decout[i];   // leave decout untouched
	  }
	  if(thisenctype == nonsystematicenc) {
		 unencode(messagebitsout);
	  }
	  else if(thisenctype == systematicenc) {	
		 extractmessagefromsyst(messagebitsout,decout);
	  }
	  else {
		 cerr << "thisenctype needs to be set" << endl;
	  }
   }
   info = messagebitsout;
   return true;
}


int polarcode::assignInitialPath()
// Tal/Vardy Algorithm 6  (some for llr and non-llr)
{

   int l = inactivePathIndices.pop();
   activePath[l] = true;
   // associate arrays with path index
   for(int lambda = 0; lambda <= n; lambda++) {
	  int s = inactiveArrayIndices[lambda].pop();
	  pathIndexToArrayIndex[lambda][l] = s;
	  arrayReferenceCount[lambda][s] = 1;
   }
   return l;
}

int polarcode::clonePath(int l)
// Tal/Vardy Algorithm 7
{
   int lp = inactivePathIndices.pop();
   activePath[lp] = true;
   // Make lp reference same arrays as l
   for(int lambda = 0; lambda <= n; lambda++) {
	  int s = pathIndexToArrayIndex[lambda][l];
	  pathIndexToArrayIndex[lambda][lp] = s;
	  arrayReferenceCount[lambda][s]++;
   }
   return lp;
}

int polarcode::clonePathllr(int l)
// Tal/Vardy Algorithm 7
{
   int lp = inactivePathIndices.pop();
   activePath[lp] = true;
   pathmetric[lp] = pathmetric[l];
   // Make lp reference same arrays as l
   for(int lambda = 0; lambda <= n; lambda++) {
	  int s = pathIndexToArrayIndex[lambda][l];
	  pathIndexToArrayIndex[lambda][lp] = s;
	  arrayReferenceCount[lambda][s]++;
   }
   return lp;
}


void polarcode::killPath(int l)
// Tal/Vardy Algorithm 8
{
   // Mark the path index l as inactive
   activePath[l] = false;
   inactivePathIndices.push(l);
   // Disassociate arrays with path index
   for(int lambda = 0; lambda <= n; lambda++) {
	  int s = pathIndexToArrayIndex[lambda][l];
	  arrayReferenceCount[lambda][s]--;
	  if(arrayReferenceCount[lambda][s] == 0) {
		 inactiveArrayIndices[lambda].push(s);
	  }
   }
}

double ** polarcode::getArrayPointer_P(int lambda, int l)
// Tal/Vardy Algorithm 9
{

   int sp;
   int s = pathIndexToArrayIndex[lambda][l];
   if(arrayReferenceCount[lambda][s] == 1) {
	  sp = s;
   }
   else {
	  sp = inactiveArrayIndices[lambda].pop();
	  double **Pto = arrayPointer_P[lambda][sp];
	  double **Pfrom = arrayPointer_P[lambda][s];
	  BITTYPE **Cto = arrayPointer_C[lambda][sp];
	  BITTYPE **Cfrom = arrayPointer_C[lambda][s];
	  
	  for(int i = 0; i < (1 << (n-lambda)); i++) {
		 Pto[i][0] = Pfrom[i][0];
		 Pto[i][1] = Pfrom[i][1];
		 Cto[i][0] = Cfrom[i][0];
		 Cto[i][1] = Cfrom[i][1];
	  }
	  arrayReferenceCount[lambda][s]--;
	  arrayReferenceCount[lambda][sp] = 1;
	  pathIndexToArrayIndex[lambda][l] = sp;
   }
   return arrayPointer_P[lambda][sp];
}


double * polarcode::getArrayPointer_llr(int lambda, int l)
// Tal/Vardy Algorithm 9
{

   int sp;
   int s = pathIndexToArrayIndex[lambda][l];
   if(arrayReferenceCount[lambda][s] == 1) {
	  sp = s;
   }
   else {
	  sp = inactiveArrayIndices[lambda].pop();
	  double *Pto = arrayPointer_llr[lambda][sp];
	  double *Pfrom = arrayPointer_llr[lambda][s];
	  BITTYPE **Cto = arrayPointer_C[lambda][sp];
	  BITTYPE **Cfrom = arrayPointer_C[lambda][s];
	  
	  for(int i = 0; i < (1 << (n-lambda)); i++) {
		 Pto[i] = Pfrom[i];
		 Cto[i][0] = Cfrom[i][0];
		 Cto[i][1] = Cfrom[i][1];
	  }
	  arrayReferenceCount[lambda][s]--;
	  arrayReferenceCount[lambda][sp] = 1;
	  pathIndexToArrayIndex[lambda][l] = sp;
   }
   return arrayPointer_llr[lambda][sp];
}


BITTYPE ** polarcode::getArrayPointer_C(int lambda, int l)
// Tal/Vardy Algorithm 9a
{
   int sp;

   int s = pathIndexToArrayIndex[lambda][l];
   if(arrayReferenceCount[lambda][s] == 1) {
	  sp = s;
   }
   else {
 	  sp = inactiveArrayIndices[lambda].pop();
	  double **Pto = arrayPointer_P[lambda][sp];
	  double **Pfrom = arrayPointer_P[lambda][s];
	  BITTYPE **Cto = arrayPointer_C[lambda][sp];
	  BITTYPE **Cfrom = arrayPointer_C[lambda][s];
	  for(int i = 0; i < (1 << (n - lambda)); i++) {
		 Pto[i][0] = Pfrom[i][0];
		 Pto[i][1] = Pfrom[i][1];
		 Cto[i][0] = Cfrom[i][0];
		 Cto[i][1] = Cfrom[i][1];
	  }
	  arrayReferenceCount[lambda][s]--;
	  arrayReferenceCount[lambda][sp] = 1;
	  pathIndexToArrayIndex[lambda][l] = sp;
   }
   return arrayPointer_C[lambda][sp];
}


BITTYPE ** polarcode::getArrayPointer_Cllr(int lambda, int l)
// Tal/Vardy Algorithm 9a
{
   int sp;

   int s = pathIndexToArrayIndex[lambda][l];
   if(arrayReferenceCount[lambda][s] == 1) {
	  sp = s;
   }
   else {
 	  sp = inactiveArrayIndices[lambda].pop();
	  double *Pto = arrayPointer_llr[lambda][sp];
	  double *Pfrom = arrayPointer_llr[lambda][s];
	  BITTYPE **Cto = arrayPointer_C[lambda][sp];
	  BITTYPE **Cfrom = arrayPointer_C[lambda][s];
	  for(int i = 0; i < (1 << (n - lambda)); i++) {
		 Pto[i] = Pfrom[i];
		 Cto[i][0] = Cfrom[i][0];
		 Cto[i][1] = Cfrom[i][1];
	  }
	  arrayReferenceCount[lambda][s]--;
	  arrayReferenceCount[lambda][sp] = 1;
	  pathIndexToArrayIndex[lambda][l] = sp;
   }
   return arrayPointer_C[lambda][sp];
}

void polarcode::recursivelyCalcP(int lambda, int phi)
// Tal/Vardy Algorithm 10
{
   if(lambda == 0) return;
   int psi = phi >> 1;
   // Recurse first if needed
   if((phi & 1) == 0) {
	  recursivelyCalcP(lambda-1, psi);
   }
   // Perform the calculation
   double sigma = 0;
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == false) continue;
	  double** Plambda = getArrayPointer_P(lambda,l);
	  double** Plambdam1 = getArrayPointer_P(lambda-1,l);
	  BITTYPE** Clambda = getArrayPointer_C(lambda,l);
	  for(int beta = 0; beta < (1 << (n - lambda)); beta++) {
		 if((phi & 1) == 0) {  // apply Equation (4) of Tal/Vardy
			for(int up = 0; up < 2; up++) {
			   double sum = 0;
			   for(int upp = 0; upp < 2; upp++) {
				  sum += Plambdam1[2*beta][up ^ upp] * Plambdam1[2*beta+1][upp];
			   } // for upp
			   Plambda[beta][up] = 0.5*sum;
			   sigma = MAX(sigma,Plambda[beta][up]);
			}  // for up
		 }
		 else { // apply Equation (5) of Tal/Vardy
			int up = Clambda[beta][0];
			for(int upp = 0; upp < 2; upp++) {
			   Plambda[beta][upp] = 0.5*Plambdam1[2*beta][up ^ upp]*
				  Plambdam1[2*beta+1][upp];
			   sigma = MAX(sigma,Plambda[beta][upp]);
			}
		 }
	  } // for beta
   } // for l

   // Normalize the probabilities
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == false) continue;
	  double **Plambda = getArrayPointer_P(lambda,l);
	  for(int beta = 0; beta < (1 << (n - lambda)); beta++) {
		 for(int u = 0; u < 2; u++) {
			Plambda[beta][u] = Plambda[beta][u]/sigma;
		 } // for u
	  } // for beta
   } // for l

}

void polarcode::recursivelyCalcllr(int lambda, int phi)
// Tal/Vardy Algorithm 10
{
   if(lambda == 0) return;
   int psi = phi >> 1;
   // Recurse first if needed
   if((phi & 1) == 0) {
	  recursivelyCalcllr(lambda-1, psi);
   }
   // Perform the calculation
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == false) continue;
	  double* LLRlambda = getArrayPointer_llr(lambda,l);
	  double* LLRlambdam1 = getArrayPointer_llr(lambda-1,l);
	  BITTYPE** Clambda = getArrayPointer_Cllr(lambda,l);
	  for(int beta = 0; beta < (1 << (n - lambda)); beta++) {
		 if((phi & 1) == 0) {  // apply Equation (4) of Tal/Vardy
			LLRlambda[beta] =
			   logdomain_sum(LLRlambdam1[2*beta]+LLRlambdam1[2*beta+1],0) - 
			   logdomain_sum(LLRlambdam1[2*beta],LLRlambdam1[2*beta+1]);
		 }
		 else { // apply Equation (5) of Tal/Vardy
			int up = Clambda[beta][0];
			LLRlambda[beta] = (1-2*up)*LLRlambdam1[2*beta] +
										LLRlambdam1[2*beta+1];
		 }
	  } // for beta
   } // for l
}


void polarcode::recursivelyUpdateC(int lambda, int phi)
// Tal/Vardy Algorithm 11
{
   int psi = phi >> 1;
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == false) continue;
	  BITTYPE ** Clambda = getArrayPointer_C(lambda,l);
	  BITTYPE ** Clambdam1 = getArrayPointer_C(lambda-1,l);
	  for(int beta = 0; beta < (1 << (n-lambda)); beta++) {
		 Clambdam1[2*beta][psi & 1] = Clambda[beta][0] ^ Clambda[beta][1];
		 Clambdam1[2*beta+1][psi & 1] = Clambda[beta][1];
	  }
   }
   if(psi & 1) {
	  recursivelyUpdateC(lambda-1,psi);
   }
}

void polarcode::recursivelyUpdateCllr(int lambda, int phi)
// Tal/Vardy Algorithm 11
{
   int psi = phi >> 1;
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == false) continue;
	  BITTYPE ** Clambda = getArrayPointer_Cllr(lambda,l);
	  BITTYPE ** Clambdam1 = getArrayPointer_Cllr(lambda-1,l);
	  for(int beta = 0; beta < (1 << (n-lambda)); beta++) {
		 Clambdam1[2*beta][psi & 1] = Clambda[beta][0] ^ Clambda[beta][1];
		 Clambdam1[2*beta+1][psi & 1] = Clambda[beta][1];
	  }
   }
   if(psi & 1) {
	  recursivelyUpdateCllr(lambda-1,psi);
   }
}


void polarcode::continuePaths_UnFrozenBit(int phi)
// Tal/Vardy Algorithm 13
{
   void sort2double_int(int num_elements, double *array, int *array2);
   int i = 0;

   // populate probForks
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == true) {
		 double **Pm = getArrayPointer_P(n,l);
		 probForks[l][0] = Pm[0][0];
		 probForks[l][1] = Pm[0][1];
		 i++;
	  }
	  else {
		 probForks[l][0] = -1;
		 probForks[l][1] = -1;
	  }
   } // for l

   int rho = MIN(2*i,L);

   // populate contForks such that contForks[l][b]is true iff
   //  probForks[l][b] is one of the rho largest entries in probForks
   //  (ties broken arbitrarily)

   memset(contForks[0],0,sizeof(bool)*2*L);
   for(int i = 0; i < 2*L; i++) {
	  idx[i] = i;
   }
   sort2double_int(2*L, probForks[0], idx);
   for(int i = 0; i < rho; i++) {
	  int j = idx[2*L - i-1];
	  int row = j / 2;
	  int col = j % 2;
	  contForks[row][col] = true;
   }

   // Kill off non-continuing paths
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == false) continue;
	  if((contForks[l][0] == false) && (contForks[l][1] == false)) {
		 killPath(l);
	  }
   } // for l
   // Then, continue relevant paths and duplicate if necessary
   for(int l = 0; l < L; l++) {
	  if((contForks[l][0] == false) && (contForks[l][1] == false)) continue;
	  BITTYPE** Cm = getArrayPointer_C(n,l);
	  if((contForks[l][0] == true) && (contForks[l][1] == true)) {
		 // both forks are good
		 Cm[0][phi & 1] = 0;
		 int lp = clonePath(l);


		 Cm = getArrayPointer_C(n,lp);
		 Cm[0][phi & 1] = 1;
		 if(do_message) {
			// copy the _message data
			for(int i = 0; i < phi; i++) {
			   arrayPointer_message[lp][i] = arrayPointer_message[l][i];
			}
			arrayPointer_message[l][phi] = 0;
			arrayPointer_message[lp][phi] = 1;
		 }
	  }
	  else {  // exactly one fork is good
		 if(contForks[l][0] == true) {
			Cm[0][phi & 1] = 0;
			if(do_message) arrayPointer_message[l][phi] = 0;
		 }
		 else {
			Cm[0][phi & 1] = 1;
			if(do_message) arrayPointer_message[l][phi] = 1;
		 }
	  }
   } // for l
}

void polarcode::continuePaths_UnFrozenBitllr(int phi)
// Tal/Vardy Algorithm 13
{
   void sort2double_int(int num_elements, double *array, int *array2);
   int i = 0;

   // populate probForks
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == true) {
		 double *LLRm = getArrayPointer_llr(n,l);
		 probForks[l][0] = -(pathmetric[l] + log(1 + exp(-LLRm[0])));
		 probForks[l][1] = -(pathmetric[l] + log(1 + exp(LLRm[0])));
		 i++;
	  }
	  else {
		 probForks[l][0] = -std::numeric_limits<double>::max();
		 probForks[l][1] = -std::numeric_limits<double>::max();
	  }
   } // for l

   int rho = MIN(2*i,L);

   // populate contForks such that contForks[l][b]is true iff
   //  probForks[l][b] is one of the rho largest entries in probForks
   //  (ties broken arbitrarily)

   memset(contForks[0],0,sizeof(bool)*2*L);
   for(int i = 0; i < 2*L; i++) {
	  idx[i] = i;
   }
   sort2double_int(2*L, probForks[0], idx);
   
   for(int i = 0; i < rho; i++) {
	  int j = idx[2*L - i-1];
	  int row = j / 2;
	  int col = j % 2;
	  contForks[row][col] = true;
   }

   // Kill off non-continuing paths
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == false) continue;
	  if((contForks[l][0] == false) && (contForks[l][1] == false)) {
		 killPath(l);
	  }
   } // for l
   // Then, continue relevant paths and duplicate if necessary
   for(int l = 0; l < L; l++) {
	  if((contForks[l][0] == false) && (contForks[l][1] == false)) continue;
	  BITTYPE** Cm = getArrayPointer_Cllr(n,l);
	  if((contForks[l][0] == true) && (contForks[l][1] == true)) {
		 // both forks are good
		 Cm[0][phi & 1] = 0;
		 int lp = clonePathllr(l);
		 Cm = getArrayPointer_Cllr(n,lp);
		 Cm[0][phi & 1] = 1;

		 double *llr_p = getArrayPointer_llr(n,l);
		 pathmetric[l] += log(1 + exp(-llr_p[0]));
		 llr_p = getArrayPointer_llr(n,lp);
		 pathmetric[lp] += log(1 + exp(llr_p[0]));

		 // copy the message data
		 if(do_message) {
			for(int i = 0; i < phi; i++) {
			   arrayPointer_message[lp][i] = arrayPointer_message[l][i];
			}
			arrayPointer_message[l][phi] = 0;
			arrayPointer_message[lp][phi] = 1;
		 }
	  }
	  else {  // exactly one fork is good
		 if(contForks[l][0] == true) {
			Cm[0][phi & 1] = 0;
			double *llr_p = getArrayPointer_llr(n,l);
			pathmetric[l] += log(1 + exp(-llr_p[0]));
			if(do_message) arrayPointer_message[l][phi] = 0;
		 }
		 else {
			Cm[0][phi & 1] = 1;
			double *llr_p = getArrayPointer_llr(n,l);
			pathmetric[l] += log(1 + exp(llr_p[0]));
			if(do_message) arrayPointer_message[l][phi] = 1;
		 }
	  }
   } // for l
}




int polarcode::findMostProbablePath(bool check_crc, encodertype thisenctype)
// in Tal/Vardy IT paper, this is absorbed into Alg. 12
{
   int lp = 0;
   double pp = 0;
   BITTYPE ** Cm;
   double ** Pm;
   bool path_with_crc_pass = false;
   unsigned int check;
   
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == false) continue;
	  
	  if(check_crc) {
		 if(thisenctype == nonsystematicenc) {
			extractmessage(messagebitsout,arrayPointer_message[l]);
			check = CRCcheck(messagebitsout,K);
		 }
		 else if(thisenctype == systematicenc) {
			// pull the message (including CRC) out of the codeword
            //  and check its parity
			BITTYPE **C0 = getArrayPointer_C(0,l);
			for(int i = 0, j=0; i < N; i++) {
			   if(polarcodedesign[i] == -1) { // unfrozen
				  messagebitsout[j++] = C0[bitreverser[i]][0];
			   }
			}
			check = CRCcheck(messagebitsout,K);
		 }
		 else {
			cerr << "thisenctype must be either systematicenc"
			   " or nonsystematicenc" << endl;
		 }
		 if(check) continue;
	  }

	  // this path checked crc
	  path_with_crc_pass = true;
	  
	  Cm = getArrayPointer_C(n,l);
	  Pm = getArrayPointer_P(n,l);
	  if(Pm[0][Cm[0][1]] > pp) {
		 lp = l;
		 pp = Pm[0][Cm[0][1]];
	  }
   }
   if(path_with_crc_pass) {		// a path with crc checked has been found
	  return lp;
   }
   else {						// just find the best w.o. regard to crc
	  return findMostProbablePath(false,thisenctype);
   }
}


int polarcode::findMostProbablePathllr(bool check_crc, encodertype thisenctype)
// in Tal/Vardy IT paper, this is absorbed into Alg. 12,
// modified for LLR computations
{
   int lp = 0;
   double pllr = std::numeric_limits<double>::max();
   bool path_with_crc_pass = false;
   unsigned int check;
   
   for(int l = 0; l < L; l++) {
	  if(activePath[l] == false) continue;

	  if(check_crc) {
		 if(thisenctype == nonsystematicenc) {
			extractmessage(messagebitsout,arrayPointer_message[l]);
			check = CRCcheck(messagebitsout,K);
		 }
		 else if(thisenctype == systematicenc) {
			// pull the message (including CRC) out of the codeword
            //  and check its parity
			BITTYPE **C0 = getArrayPointer_C(0,l);
			for(int i = 0, j=0; i < N; i++) {
			   if(polarcodedesign[i] == -1) { // unfrozen
				  messagebitsout[j++] = C0[bitreverser[i]][0];
			   }
			}
			check = CRCcheck(messagebitsout,K);
		 }
		 else {
			cerr << "thisenctype must be either systematicenc"
			   " or nonsystematicenc" << endl;
		 }
		 if(check) continue;
	  }
	  path_with_crc_pass = true;
	  
	  if(pathmetric[l] < pllr) {
		 pllr = pathmetric[l];
		 lp = l;
	  }
   }
   if(path_with_crc_pass) { // a path with crc checked has been found
	  return lp;
   }
   else {
	  return findMostProbablePathllr(false, thisenctype);
   }
}


BITTYPE *polarcode::encodesystwithCRC(BITTYPE *u)
// Encode systematically with CRC 
{
   CRCenc(u,K);

   // fill column 0 of Bsys with 0 for indices not in A
   // fill column n of Bsys with codebits for indices in A
   memset(Bsys[0],0,sizeof(BITTYPE)*N*(n+1));
   for(int i = 0, j = 0; i < N; i++) {
	  if(polarcodedesign[i] == -1) { // not frozen
		 Bsys[n][i] = u[j++];
	  }
	  else { // frozen
		 Bsys[0][i] = 0;
	  }
   }
   int s, delta, t, l, bl, kappa;
   for(int i = N-1; i >= 0; i--) {
	  if(polarcodedesign[i] == -1) { // not frozen
		 s = n-1;
		 delta = -1;
	  }
	  else { 					// frozen
		 s = 1;
		 delta = 1;
	  }
	  t = s;
	  for(int j = 0; j < n; j++) { // for each level needing computation
		 l = MIN(t, t-delta);
		 bl = (i & (1 << (n-l-1)));
		 if(bl == 0) { 			// a line with an adder at this level
			kappa = (1 << (n-1-l));
			Bsys[t][i] = Bsys[t - delta][i] ^ Bsys[t - delta][i + kappa];
		 }
		 else { 				//  a line with no adder at this level -- copy
			Bsys[t][i] = Bsys[t - delta][i];
		 }
		 t += delta;   // t = s + j*delta
	  } // for j
	  // uwork[i] = Bsys[n][i];  // copy without bit reverse permutation
	  uwork[bitreverser[i]] = Bsys[n][i];  // copy with bit reverse permutation
   } // for i

   return uwork;
}


void polarcode::setchannelfact(double cf)
{
   channelfact = cf;
}

void polarcode::computeLLR(double *y, double *llrout)
{
   for(int i = 0; i < N; i++) {
	  llrout[i] = channelfact*y[i];
   }
}


void polarcode::setAWGNchannelparams(double EbN0dB)
{
   R = double(K)/double(N);
   Ec = 1;				// coded BPSK energy
   Ecsqrt = sqrt(Ec);
   Eb = Ec/R;			// energy per bit
   EbN0 = pow(10.0, EbN0dB/10);
   // double N0, sigma2, sigma;
   N0 = Eb/(EbN0); sigma2 = N0/2;  sigma = sqrt(sigma2);
   setchannelfact(-4*sqrt(Ec)/N0);
}


double polarcode::W(double y, BITTYPE bit)
// Compute the basic likelihood function
{
   double exp(double);
   double s = Ecsqrt*(2*bit-1);
   double d = y - s;
   return exp(-(d*d)/(2*sigma2));
}


// Stuff used for CRC encode/decode

BITTYPE* polarcode::encodewithCRC(BITTYPE* u)
// Encode nonsystematically with CRC 
{

   CRCenc(u,K);
   // Merge u in with the frozen data and bit-reverse permute
   for(int i = 0, j=0; i < N; i++) {
	  if(polarcodedesign[i] == -1) { // not a frozen bit
		 uwork[bitreverser[i]] = u[j++];
	  }
	  else {					// a frozen bit
		 uwork[bitreverser[i]]= polarcodedesign[i];
	  }
   }
   for(int i = 0; i < N; i++) { // save the input message (for debug/test)
	  uworksave[i] = uwork[bitreverser[i]];
   }
   vectFn(uwork);  
   return uwork;
}

unsigned int polarcode::CRCenc(BITTYPE *data, int len)
// CRC encode the data
// len is length of data INCLUDING the 16 bits of CRC
{
   // set the last 16 bits to 0
   for(int i = len - crcsize; i < len; i++) {
	  data[i] = 0;  				// set 16 bits to 0
   }
   unsigned int state = CRCpolydiv(data,len);

   // set the last 16 bits to the parity values
   for(int i = 0; i < 16; i++) {
	  if(state & (1 << (15-i))) data[len-16+i] = 1;
	  else data[len-16+i] = 0;
   }
   return state;
}

unsigned int polarcode::CRCcheck(BITTYPE *data, int len)
// check the CRC encoding
// len is length of data INCLUDING the 16 bits of CRC
{
   return CRCpolydiv(data,len);
}

unsigned int polarcode::CRCpolydiv(BITTYPE *data,int len)
{
   unsigned int state = 0;
   unsigned char out;

   // do the polynomial division
   //   part 1: initial shift in
   int i;
   int crcn = len - 1;
   for(i = 0; i < MIN(crcp,crcn); i++) {
	  state = (state<<1) ^ data[i];
   }
   //    part 2:shift the rest in
   for( ; i <= crcn; i++) {
	  out = state >> (crcp-1);
	  state = (((state<<1)^data[i])^(crcg*out))&crcmask;
   }


   // return the parity bytes
   return state;
}


void polarcode::set_do_message()
{
   do_message = true;
   if(arrayPointer_message == 0) {
	  CALLOCMATRIX(arrayPointer_message, BITTYPE,L,N);
   }
}


void polarcode::extractmessage(BITTYPE *messageout, BITTYPE *messagein)
// take an N-bit message vector and extract the information bits, excluding
// CRC bits if any, and place in the correct order
// The message is placed in messageout, a vector of length K - crcsize
{
   for(int i = 0, j=0; i < N; i++) {
	  if(polarcodedesign[i] == -1) { // unfrozen bit
		 messageout[j++] = messagein[i];
	  }
   }
}

void polarcode::extractmessagefromsyst(BITTYPE *message, BITTYPE *codeword)
// extract the message bits from a systematically encoded codeword
// Take an N-bit systematically encoded codeword and
// extract the bits, placing the results in messageout.
{
   for(int i = 0, j= 0; i < N; i++) {
	  if(polarcodedesign[i] == -1) {  // unfrozen
		 message[j++] = codeword[bitreverser[i]];
	  }
   }
}


void polarcode::dumpdat(void)
{
   int i;
cout << "----------------" << endl;
   cout << "inactivePathIndices"; inactivePathIndices.printstack();
   printboolvec("activePath",activePath,L);
   for(i = 0; i <= n; i++) {
	  cout << "inactiveArrayIndices[" <<  i << "]: ";
	  inactiveArrayIndices[i].printstack();
	  cout << endl;
   }
   printintmat("pathIndexToArrayIndex",pathIndexToArrayIndex,n+1,L);
   printintmat("arrayReferenceCount",arrayReferenceCount,n+1,L);
cout << "----------------" << endl;
}

/*
Local Variables:
compile-command: "g++ -c polarcode.cc -std=c++11"
End:
*/


// compile-command: "g++ -c polarcode.cc -std=c++11 -stdlib=libc++"

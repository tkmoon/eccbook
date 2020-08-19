//
//
//  Program: testpolar.cc
//
//  Todd K. Moon
//  Utah State University
//
//  Date: Started Aug 15, 2018 (incorporating previous stuff)
//

#include <iostream>
#include <fstream>
#include "polarcode.h"   // conventional successive cancellation decoder
#include <math.h>
#include <cstring> // memset
#include "matalloc.h"
using namespace std;

#include "printstuff.cc"
double uran(void);
double gran(void);
void randbits(BITTYPE *u, int N);
void randnoise(double *n, double sigma, int N);
void addnoise(double *s, double *n, double *y, int N);
void bpskmodbits(BITTYPE *x, double *s, double Ecsqrt, int N);
void addscalenoise(double *s, double *n, double sigma, double *y, int N);

// data for various experiments
typedef bool (polarcode::*p2decoderfunc)(double *, BITTYPE* &, BITTYPE*&,
											 encodertype);
// p2decoderfunc is a pointer to a function (one of the decoder functions)
typedef struct experiment { 
   int K;						// actual number of message bits
   BITTYPE *message;     // message in
   BITTYPE *codeword;    // codeword
   double *s;			 // modulated signal
   polarcode *polar;     // polarcode object used for decoding
   p2decoderfunc decoderfunc;	// polarcode decoder function
   encodertype thisenctype;		// argument to encoder function
   string description;			// description of experiment
} experimenttype;
void printresults(ostream& os, unsigned int niter, experimenttype* experiments,
				  int nexp,
			 double *EbN0dBList, int neb, unsigned int** worderrcount,
			 unsigned int **wordcount, unsigned int **biterrcount,
			 unsigned int **bitcount);
void printresults2(ostream& os, unsigned int niter, experimenttype* experiments,
				   int nexp,
			 double *EbN0dBlist, int neb, unsigned int** worderrcount,
			 unsigned int **wordcount, unsigned int **biterrcount,
				   unsigned int **bitcount);

void counterrs(int expctr, int ebctr, int niter, int K, int N, 
			   unsigned int **wordcount, unsigned int **bitcount,
			   unsigned int **worderrcount, unsigned int **biterrcount,
			   BITTYPE *u1, BITTYPE *u2, BITTYPE *x1, BITTYPE *decodecw,
			   bool cwret,
			   unsigned int** worderrcountnsymsave,
			   unsigned int** biterrcountnsymsave,
unsigned int& minworderrct, unsigned int& minbiterrct);


int main(int argc, char *argv[])
{
   int i,j;

   
   int N;
   int Kp;						// K' = code dim without CRC
   int K;						// actual code dim (including CRC, if any)
   int crcsize;					// number of CRC bits  (=0 or 16)
   int L;						// number in list decoding

   int doprintvalue = 1;
   setdoprintvalue(doprintvalue);

   // allocate polarcode structures, with mfemories for all the options
   // ----------- L = 1 ------------------------------------------------
   L = 1;
   N = 2048;
   crcsize = 0;
   Kp = 1024 - crcsize;
   int Kmax = Kp; 				// the largest code dimension (no CRC)

   class polarcode polarcodeL1(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);

   polarcodeL1.set_do_message();

   // ------------L = 2 -------------------------------------------------
   L = 2;  
   crcsize = 0;  // without CRC
   Kp = 1024 - crcsize;
   K = Kp;
   polarcode polarcodeL2(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL2.set_do_message();

   L = 2;
   crcsize = 16;  // with CRC
   Kp = 1024 - crcsize;
   K = Kp;
   int Kmin = K;				// smallest code dimension (has CRC)
   polarcode polarcodeL2crc(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL2crc.set_do_message();

   // ------------- L = 4 --------------------------------------------------
   L = 4;
   crcsize = 0;  // without CRC
   Kp = 1024 - crcsize;
   K = Kp;
   polarcode polarcodeL4(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL4.set_do_message();

   L = 4;
   crcsize = 16;   // with CRC
   Kp = 1024 - crcsize;
   K = Kp;
   polarcode polarcodeL4crc(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL4crc.set_do_message();

   // ------------ L = 8 -----------------------------------------------
   L = 8;
   crcsize = 0;  // without CRC
   Kp = 1024 - crcsize;
   K = Kp;
   polarcode polarcodeL8(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL8.set_do_message();

   L = 8;
   crcsize = 16;
   Kp = 1024 - crcsize;
   K = Kp;
   polarcode polarcodeL8crc(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL8crc.set_do_message();

   // ------------ L = 16 -----------------------------------------------
   L = 16;
   crcsize = 0; // without CRC
   Kp = 1024 - crcsize;
   K = Kp;
   polarcode polarcodeL16(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL16.set_do_message();

   L = 16;
   crcsize = 16;  // with CRC
   Kp = 1024 - crcsize;
   K = Kp;
   polarcode polarcodeL16crc(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL16crc.set_do_message();
   
   // ------------ L = 32 -----------------------------------------------
   L = 32;
   crcsize = 0;  // without CRC
   Kp = 1024 - crcsize;
   K = Kp;
   polarcode polarcodeL32(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL32.set_do_message();

   L = 32;
   crcsize = 16;
   Kp = 1024 - crcsize;
   K = Kp;
   polarcode polarcodeL32crc(N,Kp,systematicenc | withcrcenc | nonsystematicenc,
					   SCdec | listllrdec | listdec, L,crcsize);
   polarcodeL32crc.set_do_message();



   polarcode polarcodelist[] {polarcodeL1, polarcodeL2, polarcodeL2crc,
		 polarcodeL4, polarcodeL4crc, polarcodeL4, polarcodeL8crc,
		 polarcodeL16, polarcodeL16crc, polarcodeL32, polarcodeL32crc};
   
   polarcode& polarcode1 = polarcodelist[0];


   int Ncoders = sizeof(polarcodelist)/sizeof(polarcodelist[0]);


   polarcodelist[0].designBhatt2(0);
   polarcodelist[0].designBhatt2(2);
   for(int i = 1; i < Ncoders; i++) {
	  polarcodelist[i].setdesign(polarcodelist[0].polarcodedesign);
   }

   BITTYPE **u;  	// message bits: u[0][*] = non-crc data (longer)
                                  // u[1][*] = crc data (shorter)
   CALLOCMATRIX(u,BITTYPE,2,Kmax);
   BITTYPE **x;	  // coded bits: x[0][*] = nonsystematic encoding, no CRC
                  //             x[1][*] = nonsystmatic with CRC
                  //             x[2][*] = systematic, no CRC
                  //             x[3][*] = systematic, with CRC

   int Nencodetypes = 4;
   CALLOCMATRIX(x,BITTYPE,Nencodetypes,N);
   BITTYPE *x1;					// pointer to codeword
   BITTYPE *u1, *u2;			// pointers to message
   
   double **s;		// BPSK symbols
   double *s1;
   CALLOCMATRIX(s,double,Nencodetypes,N);
   double *n = new double[N];		// AWGN
   double *y = new double[N];		// received signal
   
   double R = double(Kmax)/double(N); // won't quibble over slight changes in rate
   double Ec = 1;				// coded BPSK energy
   double Ecsqrt = sqrt(Ec);
   double Eb = Ec/R;
   
   double N0;

   ofstream outfile("experiments.txt");
   // Build a structure describing all the experiments
   experimenttype experiments[] {
	  // L=1
	  {Kmax, u[0], x[0], s[0], &polarcodeL1, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=1 LLR nonsyst"}, // nonsyst
	  {Kmax, u[0], x[1], s[1], &polarcodeL1, &polarcode::listdecodeLLR,
	  		systematicenc, "L=1 LLR syst"},       // syst

	  // L = 2
	  {Kmax, u[0], x[0], s[0], &polarcodeL2, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=2 LLR nonsyst"},  // nonsyst
	  {Kmax, u[0], x[1], s[1], &polarcodeL2, &polarcode::listdecodeLLR,
	  		systematicenc, "L=2 LLR syst"},     // syst  ******
	  {Kmin, u[1], x[2], s[2], &polarcodeL2crc, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=2 LLR nonsyst CRC"},  // nonsyst crc
	  {Kmin, u[1], x[3], s[3], &polarcodeL2crc, &polarcode::listdecodeLLR,
	  		systematicenc, "L=2 LLR syst CRC"},  // syst crc

      // L = 4
	  {Kmax, u[0], x[0], s[0], &polarcodeL4, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=4 LLR nonsyst"},  // nonsyst
	  {Kmax, u[0], x[1], s[1], &polarcodeL4, &polarcode::listdecodeLLR,
	  		systematicenc, "L=4 LLR syst"},     // syst  ******
	  {Kmin, u[1], x[2], s[2], &polarcodeL4crc, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=4 LLR nonsyst CRC"},  // nonsyst crc
	  {Kmin, u[1], x[3], s[3], &polarcodeL4crc, &polarcode::listdecodeLLR,
	  		systematicenc, "L=4 LLR syst CRC"},  // syst crc

      // L = 8
	  {Kmax, u[0], x[0], s[0], &polarcodeL8, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=8 LLR nonsyst"},  // nonsyst
	  {Kmax, u[0], x[1], s[1], &polarcodeL8, &polarcode::listdecodeLLR,
	  		systematicenc, "L=8 LLR syst"},     // syst   ******
	  {Kmin, u[1], x[2], s[2], &polarcodeL8crc, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=8 LLR nonsyst CRC"},  // nonsyst crc
	  {Kmin, u[1], x[3], s[3], &polarcodeL8crc, &polarcode::listdecodeLLR,
	  		systematicenc, "L=8 LLR syst CRC"},  // syst crc

      // L = 16
	  {Kmax, u[0], x[0], s[0], &polarcodeL16, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=16 LLR nonsyst"},  // nonsyst
	  {Kmax, u[0], x[1], s[1], &polarcodeL16, &polarcode::listdecodeLLR,
	  		systematicenc, "L=16 LLR syst"},     // syst ******
	  {Kmin, u[1], x[2], s[2], &polarcodeL16crc, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=16 LLR nonsyst CRC"},  // nonsyst crc
	  {Kmin, u[1], x[3], s[3], &polarcodeL16crc, &polarcode::listdecodeLLR,
	  		systematicenc, "L=16 LLR syst CRC"},  // syst crc

      // L = 32
	  {Kmax, u[0], x[0], s[0], &polarcodeL32, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=32 LLR nonsyst"},  // nonsyst
	  {Kmax, u[0], x[1], s[1], &polarcodeL32, &polarcode::listdecodeLLR,
	  		systematicenc, "L=32 LLR syst"},     // syst  *******
	  {Kmin, u[1], x[2], s[2], &polarcodeL32crc, &polarcode::listdecodeLLR,
	  		nonsystematicenc, "L=32 LLR nonsyst CRC"},  // nonsyst crc
	  {Kmin, u[1], x[3], s[3], &polarcodeL32crc, &polarcode::listdecodeLLR,
	  		systematicenc, "L=32 LLR syst CRC"},  // syst crc
		 
   };
   int nexp = sizeof(experiments)/sizeof(experiments[0]);
   cout << "number of experiments=" << nexp << endl;

   //   double EbN0dBlist[] = {2.5};	// SNRs at which to test in dB
   //   double EbN0dBlist[] = {1,2,3};
   double EbN0dBlist[] = {1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3};
   int neb = floor(sizeof(EbN0dBlist)/sizeof(double));
   double EbN0dB;
   double EbN0;
   bool *prev_decoded = new bool[neb];
   double sigma2, sigma;		// noise variance, std. dev.
   int nerrstocount = 100;
   int progressdisp = 1000;  // 10000; // 100000;
   unsigned int maxnumiter = 1000000;
   unsigned int **wordcount;   // wordcount[experiment][SNR]
   unsigned int **bitcount;   // bitcount[experiment][SNR]
   unsigned int **worderrcount;   // worderrcount[experiment][SNR]
   unsigned int **biterrcount;   // biterrcount[experiment][SNR]

   unsigned int **worderrcountnsymsave; // save number of symbols at time of update
   unsigned int **biterrcountnsymsave; // save number of symbols at time of update
   unsigned int minworderrct = maxnumiter;
   unsigned int minbiterrct = maxnumiter;
   
   CALLOCMATRIX(wordcount,unsigned int,nexp,neb);
   CALLOCMATRIX(bitcount,unsigned int,nexp,neb);
   CALLOCMATRIX(worderrcount,unsigned int,nexp,neb);
   CALLOCMATRIX(biterrcount,unsigned int,nexp,neb);
   CALLOCMATRIX(worderrcountnsymsave,unsigned int,nexp,neb);
   CALLOCMATRIX(biterrcountnsymsave,unsigned int,nexp,neb);

   //   srand(time(0));
   srand(1);

   BITTYPE *decodecw = new BITTYPE[N];;
   BITTYPE *decodemessage;
   

   bool stoploop = false;

   // set all the counters to 0
   memset(wordcount[0],0,sizeof(unsigned int)*nexp*neb);
   memset(bitcount[0],0,sizeof(unsigned int)*nexp*neb);
   memset(worderrcount[0],0,sizeof(unsigned int)*nexp*neb);
   memset(biterrcount[0],0,sizeof(unsigned int)*nexp*neb);
   memset(worderrcountnsymsave[0],0,sizeof(unsigned int)*nexp*neb);
   memset(biterrcountnsymsave[0],0,sizeof(unsigned int)*nexp*neb);
   
   unsigned int niter = 0;			// number of codewords generated
   while(!stoploop) {			// loop over information generated
	  if(niter && !(niter % progressdisp)) { // show some progress information
		 printresults2(cout, niter, experiments, nexp, EbN0dBlist, neb,
					  worderrcount, wordcount, biterrcount, bitcount);
		 printresults2(outfile, niter, experiments, nexp, EbN0dBlist, neb,
					  worderrcount, wordcount, biterrcount, bitcount);
		 cout << "minworderrct = " << minworderrct <<
			"   minbiterrct = " << minbiterrct << endl;
	  }

	  // generate message bits
	  randbits(u[0],Kmax);		// non CRC codeword
	  for(int i = 0; i < Kmin; i++) u[1][i] = u[0][i]; // CRC codeword

	  randnoise(n,1,N);			// generate noise (unscaled)

	  // encode in various ways
	  randnoise(n,1,N);		   // generate unit-variance noise
	  x1 = polarcodelist[0].encode(u[0]);    // nonsyst encode, no CRC
	  for(i = 0; i < N; i++) {x[0][i] = x1[i]; }
	  x1 = polarcodelist[0].encodesyst(u[0]); // syst encode, no CRC
	  for(i = 0; i < N; i++) {x[1][i] = x1[i]; }
	  x1 = polarcodelist[2].encodewithCRC(u[1]);  // nonsyst encode, with CRC
	  for(i = 0; i < N; i++) { x[2][i] = x1[i]; }
	  x1 = polarcodelist[2].encodesystwithCRC(u[1]); // syst encode, with CRC
	  for(i = 0; i < N; i++) { x[3][i] = x1[i]; }

	  // modulate each encoded codeword
	  for(j = 0; j < Nencodetypes; j++) {
		 bpskmodbits(x[j],s[j],Ecsqrt,N);
	  }
	  minworderrct = maxnumiter;
	  minbiterrct = maxnumiter;

	  for(int expctr = 0; expctr < nexp; expctr++) { // for each experiment

		 memset(prev_decoded,0,sizeof(bool)*neb);  // prev_decoded[*] = false;
		 
		 for(int ebctr = 0; ebctr < neb; ebctr++) {

			// first check to see if this codeword has decoded correctly
			// at lower SNR.  Assume that would also be decoded at higher SNR
			bool skipthisebno = false;
			for(int ebctr2 = 0; ebctr2 < ebctr; ebctr2++) {
			   if(prev_decoded[ebctr2]) {
				  skipthisebno = true;
				  break;
			   }
			}
			if(skipthisebno) break;   // don't do any more at this or higher SNR

			EbN0dB = EbN0dBlist[ebctr];
			EbN0 = pow(10.0, EbN0dB/10);
			N0 = Eb/(EbN0);
			sigma2 = N0/2;  sigma = sqrt(sigma2);
			experiments[expctr].polar->setAWGNchannelparams(EbN0dB);

			// run the signal through this AWGN channel to produce y
			addscalenoise(experiments[expctr].s,n,sigma,y,N);

			// call the appropriate decoder function
			bool cwret = ( (*experiments[expctr].polar).*
						 (experiments[expctr].decoderfunc))
			   (y,decodecw, u2, experiments[expctr].thisenctype);

			u1 = experiments[expctr].message;
			x1 = experiments[expctr].codeword;

			counterrs(expctr, ebctr, niter, experiments[expctr].K, N,
					  wordcount, bitcount, worderrcount, biterrcount,
					  u1, u2, x1, decodecw, cwret,
					  worderrcountnsymsave,
					  biterrcountnsymsave,
					  minworderrct, minbiterrct);

		 } // ebctr
	  } // for nexpctr

	  unsigned int minct = MIN(minworderrct,minbiterrct);
	  if(minct < maxnumiter && minct >= nerrstocount) {
		 stoploop = true;
	  }
	  niter++; // increment the number of symbols generated
	  if(niter == maxnumiter) stoploop = true;
   } // end while !stoploop

   // print everything on the way out
   printresults2(cout, niter, experiments, nexp, EbN0dBlist, neb,
				worderrcount, wordcount, biterrcount, bitcount);
   printresults2(outfile, niter, experiments, nexp, EbN0dBlist, neb,
				worderrcount, wordcount, biterrcount, bitcount);

   outfile.close();
   exit(0);
}


void printresults(ostream& os, unsigned int niter, experimenttype* experiments, int nexp,
			 double *EbN0dBlist, int neb, unsigned int** worderrcount,
			 unsigned int **wordcount, unsigned int **biterrcount,
			 unsigned int **bitcount)
{
   for(int expctr = 0; expctr < nexp; expctr++) {
	  os << "niter=" << niter << "  " << experiments[expctr].description << "  ";
	  for(int ebctr = 0; ebctr < neb; ebctr++) {
		 os << "EbN0dB=" << EbN0dBlist[ebctr] << "  worderrs=" <<
			worderrcount[expctr][ebctr] << "  words=" <<
			wordcount[expctr][ebctr] << "  worderror rate=" <<
			double(worderrcount[expctr][ebctr])/double(wordcount[expctr][ebctr])
			<< "  biterrs=" << biterrcount[expctr][ebctr] << "  bits="<<
			bitcount[expctr][ebctr] << "  biterr rate=" <<
			double(biterrcount[expctr][ebctr])/double(bitcount[expctr][ebctr])
			<< endl;
	  }
   }
}

void printresults2(ostream& os, unsigned int niter, experimenttype* experiments,
				   int nexp,
			 double *EbN0dBlist, int neb, unsigned int** worderrcount,
			 unsigned int **wordcount, unsigned int **biterrcount,
			 unsigned int **bitcount)
{
   os << "EbN0dBlist: ";
   for(int expctr = 0; expctr < nexp; expctr++) {
	  os << EbN0dBlist[expctr] << " ";
   }
   for(int expctr = 0; expctr < nexp; expctr++) {

	  os << "niter=" << niter << "  " << experiments[expctr].description << endl;
	  for(int ebctr = 0; ebctr < neb; ebctr++) {
		 os << "(" << worderrcount[expctr][ebctr] << "," << wordcount[expctr][ebctr]
			<< ")  ";
	  }
	  os << endl;

	  for(int ebctr = 0; ebctr < neb; ebctr++) {
		 os << double(worderrcount[expctr][ebctr])/double(wordcount[expctr][ebctr])
			<< "  ";
	  }
	  os << endl;
	  for(int ebctr = 0; ebctr < neb; ebctr++) {
		 os << double(biterrcount[expctr][ebctr])/double(bitcount[expctr][ebctr])
			<< "  ";
	  }
	  os << endl;
   }
}


void counterrs(int expctr, int ebctr, int niter, int K, int N,
			   unsigned int **wordcount, unsigned int **bitcount,
			   unsigned int **worderrcount, unsigned int **biterrcount,
			   BITTYPE *u1, BITTYPE *u2, BITTYPE *x1, BITTYPE *decodecw,
			   bool cwret,
			   unsigned int** worderrcountnsymsave,
			   unsigned int** biterrcountnsymsave,
			   unsigned int& minworderrct, unsigned int& minbiterrct)
{
   int i;
   
   // count how many words/bits have been generated
   wordcount[expctr][ebctr]++;
   bitcount[expctr][ebctr] += K;
   
   // check for correctness
   
   bool worderr = false;
   if(cwret) {		// if there is a decoded codeword
	  // check for word error
	  for(i = 0; i < N; i++) {
		 if(x1[i] != decodecw[i]) {
			worderrcount[expctr][ebctr]++;
			worderr = true;
			worderrcountnsymsave[expctr][ebctr] = niter;
			break;
		 }
	  }
   }
   if(worderr || !cwret) {// if worderror no nodecoded codeword
	  // count bit errors
	  int biterrsinword = 0;
	  for(i = 0; i < K; i++) {
		 if(u1[i] != u2[i]) {
			biterrsinword++;
		 }
	  }
	  if(biterrsinword) {
		 biterrcount[expctr][ebctr] += biterrsinword;
		 biterrcountnsymsave[expctr][ebctr] = niter;
	  }
	  // if !cwret && there are bit errors, increment word error count
	  if(biterrsinword && !cwret) {
		 worderrcount[expctr][ebctr]++;
		 worderr = true;
		 worderrcountnsymsave[expctr][ebctr] = niter;
	  }
   }

   if(worderrcount[expctr][ebctr])
	  minworderrct =   MIN(minworderrct, worderrcount[expctr][ebctr]);
   if(biterrcount[expctr][ebctr]) 
	  minbiterrct = MIN(minbiterrct, biterrcount[expctr][ebctr]);
}


/*
Local Variables:
compile-command: "g++ -o testpolar -g testpolar.cc polarcode.cc quicksortstuff.cc randstuff.cc -std=c++11 -lm;"
End:
*/



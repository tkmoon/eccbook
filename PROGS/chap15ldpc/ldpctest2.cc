// ldpctest2test.cc -- test the low-density parity-check code
// in AWGN channel
// Todd K. Moon

#include "ldpcdecoder.h"
#include <iostream>
#include <string>
using namespace std;
#include <math.h>
#include <csignal>
extern "C" {
#include <strings.h>  // for bzero
#include <unistd.h>
#include <time.h>
}

void sigtest(int sig);
extern int sigflag;   // set in the signal function sigtest
// This will allow you to abort computations on an SNR, in case they are
// going too long. (user beware)
// From the command line (in another terminal):  kill -30 <pid>,
// where pid is printed when the program starts running

void randvec(double *y,double a, int N, double sigma);
double gran(void);

int main()
{
   //   int maxnumloop = 1000;  // maximum number of decoding iterations
   int maxnumloop = 200;  // maximum number of decoding iterations
   string fname;
   int fileoffset;
   int printstuff = 0;
   int numerrtocount = 100;  	// count up to 100 errors
   int maxblockcount = 1000000;	// maximum number of blocks to do
   int minblockcount = 1000;    // do at least this many blocks
   // int minblockcount = 1;    // do at least this many blocks
   int randseed = 1;

   if(randseed) {
	  unsigned int seedtime = (unsigned int)time(NULL);
	  srand(seedtime);
   }
   
	  

   //-------------------------------------------------------------
   // Rate 1/2 data:
   // fname = "A1-2.txt"; fileoffset = 1;
   // double EbN0dBlist[] = {1, 1.1}; // {1, 1.1, 1.2, 1.3, 1.4, 1.5};  // , 1.6};

   //-------------------------------------------------------------
//tkm
   // N=1057 from MacKay website
   // fname = "1057-244-1.txt"; fileoffset = 1;
   // double EbN0dBlist[] = {3}; // {1.5, 1.8};  // {4,4.5};

   // ----------------------------------------------------------------
   // Rate 1/4 data:
   // fname = "A1-4strip.txt"; fileoffset = 1;
   // double EbN0dBlist[] = {.4,.5,.6,.7,.8,.9,1}; //,1.1,1.2,1.3,1.4};
   
   // -------------------------------------------------------------
   // Rate 1/3 data:
   // fname = "A1-3strip.txt"; fileoffset = 1;
   // double EbN0dBlist[] = {0.4,0.5,0.6,0.7,0.8,0.9,1};

   // Small test matrix data
    // fname = "Agall.txt"; fileoffset = 1;
    // double EbN0dBlist[] = {4,4.5,5,5.5,6,6.5,7,7.5,8,8.5};

   // -------------------------------------------------------------
   // Test matrix from chapter
//tkm
   // fname = "moonA.txt"; fileoffset = 0;
   // double EbN0dBlist[] = {1};

   // fname = "HRSg4r32.txt"; fileoffset = 1;
   // double EbN0dBlist[] = {3.3, 3.5, 3.7};

   // Matrix designed using PEG
   fname = "PEGN4000.txt"; fileoffset = 0;
   double EbN0dBlist[] = {1,1.2};


   // -------------------------------------------------------------
   // Hamming code
   // fname = "Ahamming74.txt";
   // fileoffset = 1;
   // double EbN0dBlist[] = {2,4,5,6,7,8};
   // printstuff = 0;


   int neb = sizeof(EbN0dBlist)/sizeof(double); // number of steps

// tkm
   cout << "reading decoder file" << endl;
  LDPCDECODER ldpcdecoder(fname,fileoffset, DOPROBDECODE); 
  //   LDPCDECODER ldpcdecoder(fname,fileoffset, DOPROBDECODE | LPDECODE ); 
  // LDPCDECODER ldpcdecoder(fname,fileoffset, DOPROBDECODE | DC1 ); 
   //   LDPCDECODER ldpcdecoder(fname,fileoffset, DOPROBDECODE | WEIGHTEDBITFLIP | MODWEIGHTEDBITFLIP |MULTIGDBITFLIP);
   // LDPCDECODER ldpcdecoder(fname,fileoffset, DOPROBDECODE | DOAMINSTARDECODE
   // 						   | DORCBPDECODE);

   // LDPCDECODER ldpcdecoder(fname,fileoffset, DOPROBDECODE | DOLOGLIKEDECODE |
   // 						   DOMINSUMDECODE | DOPHIDECODE | DOQPHIDECODE);
   // LDPCDECODER ldpcdecoder(fname,fileoffset,DOLOGLIKEDECODE|DOMINSUMDECODE);
   // LDPCDECODER ldpcdecoder(fname,fileoffset,DOPROBDECODE | DOLOGLIKEDECODE
   // 						   | DOMINSUMDECODE);

cout << "decoder file read in" << endl;

   // construct the decoder object
   int ndec = ldpcdecoder.numdecodetype;
   // set up a signal handler to break out of loops under user
   // control

   cout << "pid=" << getpid() << "    siguser1=" << SIGUSR1 << endl; 
   void (*prev_fn)(int);   // pointer to previous function
   prev_fn = signal(SIGUSR1,sigtest);  // sigtest is in ldpcdecoder.cc

   int i;
   unsigned long int maxcount; // maximum number of blocks to code
   double R = double(ldpcdecoder.K)/double(ldpcdecoder.N); // code rate
   cout << "Code rate=" << R << endl;


   double *y = new double[ldpcdecoder.N]; // vector of channel outputs
   double *p = new double[ldpcdecoder.N]; // vector of channel posterior probs

   double EbN0dB, EbN0, sigma2, sigma;
   int decval;					// return value from decoder

   // variables to hold decoder statistics
   unsigned long int ncount, *ndederrcount, *nundederrcount, nblockcount;
   unsigned long int *toterrcount;
   unsigned long int mintoterrcount; // smallest total error count
   // among all decoders
   unsigned long int *decfailcount; // count number of decoder failures
   unsigned long int *decsuccesscount; // number of decoder successes
   unsigned long *numdeciter;	// number of iterations of decoder
								// used to compute average number
   unsigned long int *maxnumdeciter; // maximum number of iterations on decode
   unsigned long int *numloops;	// number of loops to do
   numloops = new unsigned long int[ndec];
   ndederrcount = new unsigned long int[ndec];
   nundederrcount = new unsigned long int[ndec];
   toterrcount = new unsigned long int[ndec];
   decfailcount = new unsigned long int[ndec];
   decsuccesscount = new unsigned long int[ndec];
   numdeciter = new unsigned long int[ndec];
   maxnumdeciter = new unsigned long int[ndec];

   double **durations;  // durations (in seconds) of decoding
   CALLOCMATRIX(durations,double,neb,ndec);
   
   double *ebnolist = new double[neb];   // list of EbN0 values
   unsigned long int **nundedlist;       // number of undedcoded samples
   CALLOCMATRIX(nundedlist,unsigned long int,ndec,neb);
   unsigned long int **ndedlist;	// number of decoded samples
   CALLOCMATRIX(ndedlist,unsigned long int,ndec,neb);
   double **errlist;    // error probability
   CALLOCMATRIX(errlist,double,ndec,neb);
   unsigned long int **decfaillist;	     // decoder failures
   CALLOCMATRIX(decfaillist,unsigned long int, ndec,neb);
   unsigned long int **decsuccesslist;   // decoder successes
   CALLOCMATRIX(decsuccesslist,unsigned long int, ndec, neb);
   double **numdeciteravglist; // average number of decoder iterations
   CALLOCMATRIX(numdeciteravglist,double, ndec, neb);
   unsigned long int **maxnumdeciterlist;// maximum number of decoder iter's
   CALLOCMATRIX(maxnumdeciterlist,unsigned long int, ndec, neb);
   double **werlist;  // word error rate
   CALLOCMATRIX(werlist,double, ndec, neb);
   unsigned long int ncodewords;
   unsigned long int *worderrcount;
   worderrcount = new unsigned long int[ndec];
   bool allbitsinword;			// true if all decoded bits in codeword correct

   // ldpcdecoder.dosave1("savetest1.txt");
   // ldpcdecoder.dosave2("savetest2.txt");

   double a = -1;
   ldpcdecoder.setsigamp(a);			// set signal amplitude = 1 = sqrt(Ec)
   bzero(durations[0],neb*ndec*sizeof(double));
   int showprogress = 0;
   int ctr = 0;
   for(ctr=0; ctr < neb; ctr++) { // loop over various SNR values
	  if(showprogress) cout << endl;
	  EbN0dB = EbN0dBlist[ctr];	// get SNR in DB
	  EbN0 = pow(10.,EbN0dB/10.);  // get SNR
	  sigma2 = 1/(2*R*EbN0);       // compute noise variance (Ec = 1)
	  ldpcdecoder.setsigma2(sigma2);
	  sigma = sqrt(sigma2);        // and standard deviation
	  cout << "EbN0dB=" << EbN0dB << "  EbN0=" << EbN0 << "  sigma2="
		   << sigma2 << endl;
	  // initialize the variables for this SNR
	  nblockcount = 0;			// number of blocks processed
	  ncount = 0;				// total number of bits
	  ncodewords = 0;			// total number of codewords

	  bzero(ndederrcount,ndec*sizeof(unsigned long int));
	  // number of detected errors
	  bzero(nundederrcount,ndec*sizeof(unsigned long int));
	  // number of undetected errors
	  bzero(toterrcount,ndec*sizeof(unsigned long int));
	  // total number of errors
	  bzero(decfailcount,ndec*sizeof(unsigned long int));
	  // count number of decoding failures
	  bzero(decsuccesscount,ndec*sizeof(unsigned long int));
	  // count number of decoding successes
	  bzero(decsuccesscount,ndec*sizeof(unsigned long int));
	  // number of decoding iterations
	  bzero(maxnumdeciter,ndec*sizeof(unsigned long int));
	  // maximum number of decoding iterations
	  bzero(worderrcount,ndec*sizeof(unsigned long int));
	  // word error count
	  mintoterrcount = 0;

	  while(1) {
		 // count until a specified # of errors, and not too many blocks
		 nblockcount++;
		 if(showprogress) {
			if(!(nblockcount % 100)) cout << nblockcount << " ";
			else cout << "." << flush;
		 }
		 randvec(y,a,ldpcdecoder.N,sigma);  // Set up received vector, 
		 // assume all 0 codeword transmitted
		 
		 // Data for the example of the chapter
//tkm
		 // ldpcdecoder.setsigamp(-2);	// set signal amplitude
		 // ldpcdecoder.setsigma2(2);  // set noise variance
		 // y[0] = -.63; y[1] = -0.83;  y[2] = -0.73;  y[3] = -0.04;
		 // y[4] = 0.1; // error sign
		 // //y[4] = -0.1;  // correct sign
		 // y[5] = 0.95;  y[6] = -0.76;  y[7] = 0.66; y[8] = -0.55; y[9] = 0.58;
		 // minblockcount = 0; numerrtocount = 0;  maxblockcount = 1;
		 // printstuff = 0;
		 // test the bitflip:
		 // y[0] = -1; y[1] = 1; y[2] = -1; y[3] = 1; y[4] = -1;
		 // y[5] = 1; y[6] = -1; y[7] = 1; y[8] = -1; y[9] = 1;

		 bzero(numloops,ndec*sizeof(unsigned long int));
		 ncodewords++;				  // increment the number of codewords
		 ncount += ldpcdecoder.N;    // accumulate the total number of bits

		 
		 // decode loops over all the different kinds of decoders
		 decval = ldpcdecoder.decode(y , printstuff, maxnumloop, numloops,
									 durations[ctr]);

		 
		 // devcal=1 in decoder type bit position if decoded correctly
		 for(int ndecoders = 0; ndecoders < ndec; ndecoders++){
			// loop over all decoder types

			if(!(ldpcdecoder.decodetypelist[ndecoders] & decval)) {
			   // did not decode for this decoder successfully:
			   // count detected errors
			   decfailcount[ndecoders]++; // count number of decoder failures
			   // count the detected bits in error
			   allbitsinword = true;
			   for(int i = 0; i < ldpcdecoder.N; i++) {
				  if(ldpcdecoder.c[ndecoders][i] != 0) {
					 allbitsinword = false;
					 ndederrcount[ndecoders]++;
				  }
			   }
			}
			else { // did decode for this decoder successfully:
			   decsuccesscount[ndecoders]++; //count number of decoder successes
			   // accumulate number of decoder iterations
			   numdeciter[ndecoders] += numloops[ndecoders]; 
			   // save maximum number of decoder iterations
			   if(numloops[ndecoders] > maxnumdeciter[ndecoders])
				  maxnumdeciter[ndecoders] = numloops[ndecoders];
			   // count undetected errors
			   allbitsinword = true;
			   for(int i = 0; i < ldpcdecoder.N; i++) {
				  if(ldpcdecoder.c[ndecoders][i] != 0) {
					 allbitsinword = false;
					 nundederrcount[ndecoders]++;
				  }
			   }
			}
			toterrcount[ndecoders] = ndederrcount[ndecoders] +
			   nundederrcount[ndecoders]; // total number of bit errors
			if(!allbitsinword) worderrcount[ndecoders]++;
		 }  // for ndecoders

		 // find the decoder with the minimum number of errors
		 mintoterrcount = toterrcount[0];
		 for(int ndecoders = 1; ndecoders < ndec; ndecoders++) {
			if(toterrcount[ndecoders] < mintoterrcount) {
			   mintoterrcount = toterrcount[ndecoders];
			}
		 }			

		 if(nblockcount >= minblockcount) {  // contemplate breaking
			if(mintoterrcount > numerrtocount || (nblockcount > maxblockcount))
			   break;
		 }
		 if(sigflag) {
			cout << "Breaking out of this SNR" << endl;
			sigflag = 0;
			break;
		 }
		 if(ldpcdecoder.savestuff1 || ldpcdecoder.savestuff1) {
			exit(0);
		 }
	  } // end while(1)
	  cout << "nblockcount=" << nblockcount << "  total bits="  <<
		 ncount << endl;
	  for(int ndecoders = 0; ndecoders < ndec; ndecoders++) {
cout << "yo0" << endl;
		 durations[ctr][ndecoders] /= nblockcount;
		 cout << "Decoder type: " <<
			// 			ldpcdecoder.decodenames[ndecoders] << endl;
			decodernameslist[ldpcdecoder.decodetypenum[ndecoders]]
			  << endl;

cout << "yo1" << endl;
cout << "ncount=" << ncount << endl;
cout << "ncount+1=" << ncount + 1 << endl;
 long unsigned int junk = ncount;
 cout << "junk=" << junk;
 double djunk = double(junk);
 cout << "djunk=" << djunk << endl;
 
// ncount = 5;
// double dncount = ncount;
// cout << "dncount=" << dncount << endl;
// cout << "ncount=" << ncount << "   double(ncount)=" << double(ncount) << endl;
cout << "ndec=" << ndec << "  ndecoders=" << ndecoders << endl;
cout << "ndederrcount[ndecoders]=" << ndederrcount[ndecoders] << endl;
cout << "nundederrcount[ndecoders]=" << nundederrcount[ndecoders] << endl;
cout << "toterrcount[ndecoders]=" << toterrcount[ndecoders] << endl;

cout << "errorrate=" << double(toterrcount[ndecoders])/double(ncount) << endl;
cout << "worderrorrate=" << double(worderrcount[ndecoders])/double(ncodewords) << endl;
          cout << "#detected errors=" << 
			ndederrcount[ndecoders]
			  << "  # undet. errors=" << nundederrcount[ndecoders] << 
			"  toterrcount=" << toterrcount[ndecoders] << 
			"  errorrate=" <<
			double(toterrcount[ndecoders])/double(ncount) << endl;
cout << "yo2" << endl;
		 cout << "worderrrate=" << double(worderrcount[ndecoders])/
			double(ncodewords) << endl;
		 if(decsuccesscount[ndecoders])
			cout << "avgdeciter=" <<
			   double(numdeciter[ndecoders])/double(decsuccesscount[ndecoders]);
cout << "yo3" << endl;
		 cout <<"   maxdeciter=" << maxnumdeciter[ndecoders];
		 cout << " decsuccesscount=" << decsuccesscount[ndecoders] << endl;
		 cout << " durations=" << durations[ctr][ndecoders] << endl;

cout << "yo4" << endl;
		 ebnolist[ctr] = EbN0;
		 nundedlist[ndecoders][ctr] = nundederrcount[ndecoders];
		 ndedlist[ndecoders][ctr] = ndederrcount[ndecoders];
cout << "yo5" << endl;
		 errlist[ndecoders][ctr] =
			(double)(toterrcount[ndecoders])/(double)ncount;
cout << "yo6" << endl;
		 decfaillist[ndecoders][ctr] = decfailcount[ndecoders];
		 decsuccesslist[ndecoders][ctr] = decsuccesscount[ndecoders];
		 werlist[ndecoders][ctr] = double(worderrcount[ndecoders])/
										  double(ncodewords);
cout << "yo7" << endl;
		 if(decsuccesscount[ndecoders])
			numdeciteravglist[ndecoders][ctr] =
			   (double)numdeciter[ndecoders]/(double)decsuccesscount[ndecoders];
		 else
			numdeciteravglist[ndecoders][ctr] = 0;
		 maxnumdeciterlist[ndecoders][ctr] = maxnumdeciter[ndecoders];
cout << "yo8" << endl;
	  }	// for ndecoders
	  if(ldpcdecoder.decodetype & MULTIGDBITFLIP) {
		 printf("multibit gdbitflip: fmin=%g  fmax=%g\n",ldpcdecoder.fmin,
				ldpcdecoder.fmax);
	  }
cout << "yo1" << endl;
   } // for ctr (loop over ebn0)
   // print out summary information
   cout << "--------------------------------------" << endl;
   cout << "R=" << R << endl;
   for(int ndecoders = 0; ndecoders < ndec; ndecoders++){
	  //	  cout << ldpcdecoder.decodenames[ndecoders] << endl;
	  cout << decodernameslist[ldpcdecoder.decodetypenum[ndecoders]]
		   << endl;

	  cout<<"EbN0" << "\t" << "nunded" << "\t" << "nded" << "\t" << "err"
		  << "\t\t" << "#fail" << "\t" << "#succ" << "\t" << "wer" <<
		 "\t\t\t" << "iteravg"
		  << "\t" << "itermax" << "\t" << "dur(s)" << endl;

	  for(i = 0; i < ctr; i++) {
		 cout << EbN0dBlist[i] <<"\t" << 
			nundedlist[ndecoders][i] << "\t" // number undetected errors
			  << ndedlist[ndecoders][i]         // 
			  << " \t" << errlist[ndecoders][i] << "\t" <<
			decfaillist[ndecoders][i] << "\t" << 
			decsuccesslist[ndecoders][i] << "\t" <<
			werlist[ndecoders][i] << "\t\t" << 
			numdeciteravglist[ndecoders][i] << "\t" << 
			maxnumdeciterlist[ndecoders][i] << "\t" <<
			durations[i][ndecoders]
			  << endl;
	  }
   }

   FREEMATRIX(nundedlist);  FREEMATRIX(ndedlist); FREEMATRIX(errlist);
   FREEMATRIX(decfaillist); FREEMATRIX(decsuccesslist);
   FREEMATRIX(numdeciteravglist); FREEMATRIX(maxnumdeciterlist);
   FREEMATRIX(durations);
   FREEMATRIX(werlist);
   
}


void randvec(double *y,double a, int N, double sigma)
// fill y with mean= -1 and variance sigma^2 Gaussian
// This corresponds to the all-zero vector
{
   int i;
   for(i = 0; i < N; i++) {
	  y[i] = a + sigma*gran();
   }
}

/* If LP decoding is not used, the following will compile/link
compile-command:"g++ -o ldpctest2 -g -std=c++11 ldpctest2.cc  ldpcdecoder.cc gran.cc"
*/

/*
Local Variables:
compile-command:"g++ -o ldpctest2 -g -std=c++11 ldpctest2.cc  ldpcdecoder.cc gran.cc -lglpk"
End:
*/

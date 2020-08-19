// ldpcdec.cc -- LDPC decoder class definitions
// Todd K. Moon

#include "ldpcdec.h"
#include <fstream>
#include <iostream>
using namespace std;
#include <iomanip>
#include <stdlib.h>
#include <math.h>

const double TINYDIV = 1.e-10;
// const double EPS = 0;
// const double CLIPONE = 1-EPS;
// #define CLIPCHECK(q0,q1,one) if(q0>one){q0=one;q1=1-one;} \
//                              else if(q1>one){q1=one;q0=1-one;}

LDPCDEC::LDPCDEC(const string fname, int offset, int indecodetype)
// fname = file containing sparse description of code
// offset = index offset.  0: C-like indexing   1: matlab-like indexing
// decodetype = 1 for probability decoder, 2 for log likelihood decoder
{
   int n,m,k,d;
   const int maxline = 1024;	// maximum length of line expected 
   char line[maxline];

   ifstream infile(fname);
   if(!infile) {
	  cerr << "Error: unable to open input file" << fname << endl;
	  exit(-1);
   }
   decodetype = indecodetype;

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
}

void
LDPCDEC::allocdecodedat(void)
// allocate memory used in the decoder
{
   int i,j;

   c = new unsigned char[N];
   na = new int[N];
   if(decodetype & 1 ) {			// probability decoder
	  pn = new double[N];
	  deltaq = new double[maxrowwt];
	  CALLOCMATRIX(r1,double,maxcolwt,N);
	  CALLOCMATRIX(r0,double,maxcolwt,N);
	  CALLOCMATRIX(q1,double,maxcolwt,N);
	  CALLOCMATRIX(q0,double,maxcolwt,N);
	  CALLOCMATRIX(deltar,double,maxcolwt,N);
	  q0p = new double[N];
	  q1p = new double[N];
   }
   if(decodetype & 2) {					// loglikelhood decoder
	  rn = new double [N];
	  lambda = new double[N];
	  CALLOCMATRIX(u,double,maxcolwt,N);
	  CALLOCMATRIX(oldu,double,maxcolwt,N);
   }
}

void
LDPCDEC::freedecodedat(void)
{
   delete [] c;
   delete [] na;
   if(decodetype & 1) {
	  delete[] q1p;
	  delete[] q0p;
	  FREEMATRIX(deltar);
	  FREEMATRIX(q0);
	  FREEMATRIX(q1);
	  FREEMATRIX(r0);
	  FREEMATRIX(r1);
	  delete[] deltaq;
	  delete [] pn;
   }
   if(decodetype & 2) {
	  FREEMATRIX(oldu);
	  FREEMATRIX(u);
	  delete[] lambda;
	  delete[] rn;
   }

}

   
int
LDPCDEC::decode(double *y, int inprintstuff, int maxnumloop, int &numloops)
// pn = channel posterior probabilities (one for each N)
// printstuff - set to print intermediate results
// maxnumloop = maximum number of decoding iterations
// numloops = number of decoding iterations actually used
// if numloops is a positive integer on input, then
// the decoder does exactly that many iterations.  
// if numloops is 0 on input, then the decoder
// proceeds until Hx = 0 (parities all check) or until maxnumloop iterations
//
// returns: 1 if decoding succeeds; 0 for a decoding failure
{
   int i,l,k,m,n,row;
   double prod, prod0, prod1,alpha;
   int paritycheck = 0;
   char z;
   int loopcount = 0;
   double sum;
   int idx;

   printstuff = inprintstuff;

   if(decodetype  & 1) {
	  probdecodeinit(y);
   }
   if(decodetype & 2) {
	  logdecodeinit(y);
   }

   if(numloops > 0) maxnumloop = numloops;

   do {  // loop until completion of decoding, or loop count
	  loopcount++;
	  if(printstuff) cout << "Iteration: " << loopcount << endl;

	  if(decodetype & 1) {
		 probdecode();  // call the actual decoding function
	  }
	  if(decodetype & 2) {
		 lldecode();    // call the actual decoding function
	  }

	  // Check the parity condition 
	  if(numloops == 0) {	// if break based on parity, not on specified 
		 // number of loops
		 paritycheck = 1;
		 for(m = 0; m < M; m++) { // check all parities
			z = 0;
			for(l = 0; l < Nmlen[m]; l++) {
			   z += c[Nm[m][l]];
			}
			z %= 2;
			if(z) {
			   // Parity check fails.  We could bail out of check at this 
			   //   point.
			   paritycheck = 0;
			   break;	// break out of parity check loop
			}
		 }
		 if(paritycheck) break;	// break out of decode loop
	  }
   } while(loopcount < maxnumloop);

   // If not previously done, check parity
   if(numloops > 0) {  // in this case, we haven't checked parity
	  paritycheck = 1;
	  for(m = 0; m < M; m++) { // check all parities
		 z = 0;
		 for(l = 0; l < Nmlen[m]; l++) {
			z += c[Nm[m][l]];
		 }
		 z %= 2;
		 if(z) {
			// Parity check fails.  We could bail out of check at this 
			//   point.
			paritycheck = 0;
			break;	// break out of parity check loop
		 }
	  }
   }
   numloops = loopcount;
   return(paritycheck);
}



void LDPCDEC::probdecodeinit(double *y)
   // convert to posterior probabilities
{
   int i,m,n;
   double L;
   if(printstuff) VECDUMP(y,N);

   L = -2*a/sigma2;
   // compute channel posterior probabilities
   for(i = 0; i < N; i++) {
	  pn[i] = 1/(1+exp(L*y[i]));
   }
   if(printstuff) VECDUMP(pn,N);
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 q1[m][n] = pn[n];
	  }
   }
   if(printstuff) printsparse("q(init)",q1);
}

void LDPCDEC::logdecodeinit(double *y)
{
   int i,j;
   double Lc;
   Lc = 2*a/sigma2;
   for(i = 0; i < N; i++) {
	  rn[i] = Lc*y[i];
	  lambda[i] = rn[i];
   }
   for(i = 0; i < maxcolwt; i++) {
	  for(j = 0; j < N; j++) {
		 oldu[i][j] = 0;
	  }
   }
}


void LDPCDEC::probdecode()
{
   int i,row,l,idx,n,k;
   double prod,prod0,prod1,sum,alpha;
   double r;

   for(i = 0; i < N; i++) {
	  na[i] = 0;				// clear the "number above" offset
	  // used for the sparse representation
   }

   // Horizontal step
   for(row = 0; row < M; row++) {
	  // copy the data on this row into a temporary array
	  for(l = 0; l < Nmlen[row]; l++) {
		 idx = Nm[row][l];	  // column index of lth element on this row
		 deltaq[l] = 1-2*q1[na[idx]][idx]; // compute delta q
	  }
	  // work over nonzero elements of this row
	  for(l = 0; l < Nmlen[row]; l++) {
		 prod = 1;
		 for(k = 0; k < Nmlen[row]; k++) {
			if(k==l) continue;
			prod *= deltaq[k];
		 }
		 // assign back into sparse structure
		 idx = Nm[row][l];
		 r1[na[idx]][idx] = (1-prod)/2;
		 r0[na[idx]++][idx] = (1+prod)/2;
	  }
   }
   if(printstuff) printsparse("r1",r1);
	  
   // Vertical step
   for(n = 0; n < N; n++) {  // work over each column
	  prod0 = 1-pn[n];  prod1 = pn[n];
	  for(l = 0; l < Mnlen[n]; l++) { // compute the pseudoposteriors
		 prod0 *= r0[l][n];
		 prod1 *= r1[l][n];
	  }
	  if(prod0>prod1) {
		 r = prod1/prod0;		// <1
		 q0p[n] = 1/(1+r);
		 q1p[n] = r/(r+1);
	  }
	  else {
		 r = prod0/prod1;		// < 1
		 q0p[n] = r/(r+1);
		 q1p[n] = 1/(1+r);
	  }
	  // now get the others by dividing out the right factor
	  for(l = 0; l < Mnlen[n]; l++) {
		 prod0 = q0p[n]/r0[l][n];
		 prod1 = q1p[n]/r1[l][n];
		 if(prod0>prod1) {
			r = prod1/prod0;		// <1
			q0[l][n] = 1/(1+r);
			q1[l][n] = r/(r+1);
		 }
		 else {
			r = prod0/prod1;		// < 1
			q0[l][n] = r/(r+1);
			q1[l][n] = 1/(1+r);
		 }
	  }
	  // do the decoding
	  if(q1p[n] > 0.5) c[n] = 1; else c[n] = 0;
   }
   
   if(printstuff) printsparse("q1",q1);
   if(printstuff) VECDUMP(q1p,N);
   if(printstuff) VECDUMP2(c,N,int);
}


void LDPCDEC::lldecode(void)
{
   int m,n1,n2;
   int l,k,idx,n,i;
   double prod;

   for(i = 0; i < N; i++) {
	  na[i] = 0;				// clear the "number above" offset
	  // used for the sparse representation
   }

   // Check node update
   for(m = 0; m < M; m++) {
	  // work over nonzero elements of this row
	  for(n1 = 0; n1 < Nmlen[m]; n1++) {
		 n = Nm[m][n1];
		 prod = 1;
		 for(n2 = 0; n2 < Nmlen[m]; n2++) {
			if(n1==n2) continue;
			i = Nm[m][n2];
			prod *= tanh((-lambda[i] + oldu[na[i]][i])/2);
		 }
		 u[na[n]][n] = -2*atanh(prod);
	  }
	  for(n1 = 0; n1 < Nmlen[m]; n1++) {
		 n = Nm[m][n1];
		 na[n]++;
	  }
   }
   // copy u into oldu
   for(n = 0; n < N; n++) {
	  for(m = 0; m < Mnlen[n]; m++) {
		 oldu[m][n] = u[m][n];
	  }
   }

   if(printstuff) printsparse("u",u);
   
   // Bit node update
   for(n = 0; n < N; n++) {  // work over each column
	  lambda[n] = rn[n];
	  for(l = 0; l < Mnlen[n]; l++) {
		 lambda[n] += u[l][n];
	  }
	  if(lambda[n] >= 0) c[n] = 1; else c[n] = 0;
   }
   if(printstuff) VECDUMP(lambda,N);
   if(printstuff) VECDUMP2(c,N,int);
}

void LDPCDEC::printsparse(const string& name,double **sparsemat)
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
			// printf("%.2f\t",0.0);
		 }
		 // print the new number
		 cout << setw(5) << sparsemat[na[n]++][n] << "\t";
		 // printf("%.2f\t",sparsemat[na[n]++][n]);
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < N; n1++) {
		 cout << setw(5) << "0.0\t";
		 // printf("%.2f\t",0.0);
	  }
	  cout << endl;
	  // printf("\n");
   }
   delete[] na;
   cout.precision(sp);
}

void
LDPCDEC::printsparseA(void)  // print the sparse binary A matrix
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
			// printf("0 ");
		 }
		 // print the new number
		 cout << "1 "; 
		 // printf("1 ");
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < N; n1++) {
		 cout << "0 ";
		 // printf("0 ");
	  }
	  cout << endl;
	  // printf("\n");
   }
}


/*
Local Variables:
compile-command: "g++ -c -g ldpcdec.cc -Wno-deprecated"
End:
*/

//
//
//  Program: testpolarSCD.cc
//
//  Todd K. Moon
//  Utah State University
//
//  Date: Started June 30, 2018
//

#include <iostream>
#include "polarcode.h"

#include <math.h>
#include "matalloc.h"
using namespace std;

#include "printstuff.cc"
double uran(void);
double gran(void);
void randbits(BITTYPE *u, int N);
void randnoise(double *n, double sigma, int N);
void addnoise(double *s, double *n, double *y, int N);
void bpskmodbits(BITTYPE *x, double *s, double Ecsqrt, int N);

#define sign(x) (x >= 0 ? 1 : -1)

int main(int argc, char *argv[])
{
   int N;
   int K;
   
   N = 16;
   K = 8;
   int M = 1000000;    // number of bits to simulate for MC design

   polarcode polarcode(N,K);
   polarcode.designBhatt2(0);
   printsignedcharvec("Design Bhatt",polarcode.polarcodedesign,N);
   SIGNEDBITTYPE *savepolarBhatt = new SIGNEDBITTYPE[N];
   for(int i = 0; i < N; i++) savepolarBhatt[i] = polarcode.polarcodedesign[i];
   polarcode.designMC(0,M);
   printsignedcharvec("Design MC",polarcode.polarcodedesign,N);
   for(int i = 0; i < N; i++) {
   	  if(savepolarBhatt[i] != polarcode.polarcodedesign[i]) {
   		 cout << "designs different at i=" << i << "  Bhatt=" <<
   			int(savepolarBhatt[i]) << "  MC=" <<
   			int(polarcode.polarcodedesign[i]) << endl;
   	  }
   }
   exit(0);

}


/*
Local Variables:
compile-command: "clang++ -o testpolardesign -g testpolardesign.cc polarcode.cc quicksortstuff.cc randstuff.cc -std=c++11 -stdlib=libc++ -lm;"
End:
*/



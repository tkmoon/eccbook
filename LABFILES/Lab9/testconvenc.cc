//  Program: testconvenc.cc
//  Todd K. Moon

#include <iostream>
using namespace std;

#include "BinConvFIR.h"
#include "BinConvIIR.h"

#include "matalloc.h"
void printconvtest(BinConv &conv, int in1, int Nimpulse);

int main()
{

   // First test:
   int k1 = 1;					// one inputs
   int n1 = 2;					// two outputs
   unsigned int **g1;
   CALLOCMATRIX(g1,unsigned int, k1,n1);
   // G = [D^2+1 D^2+D+1 ]
   g1[0][0] = 5;  g1[0][1] = 7;  // first row of G
   int p1[] = {2};        // degrees of rows of G
   BinConvFIR conv1(k1,n1,p1,g1);

   int Nimpulse = 10;			// number of samples of impulse response
   int in1 = 1;					// input number to test

   cout << "First Test: " << endl;
   printconvtest(conv1,in1, Nimpulse);
   cout << endl << endl << endl;

   // Second test:
   // G = [D^2+D+1 D^2  1+D
   //      D 1   0]
   unsigned int **g2;
   int k2 = 2;					// two inputs
   int n2 = 3;					// three outputs
   CALLOCMATRIX(g2,unsigned int, k2,n2);
   g2[0][0] = 7;  g2[0][1] = 4;  g2[0][2] = 3;  // first row of G
   g2[1][0] = 2;  g2[1][1] = 1;  g2[1][2] = 0;  // second row of G
   int p2[] = {2,1};        // degrees of rows of G
   BinConvFIR conv2(k2,n2,p2,g2);
   int in2 = 2;					// input number to test
                                // set in2=1 or 2 or 3 to test impulse response
                                // of first, second, or both inputs
   cout << "Second Test: " << endl;
   printconvtest(conv2,in2, Nimpulse);
   cout << endl << endl << endl;

   // Third test: Test the IIR stuff
   int k3 = 2;					// two inputs
   int n3 = 3;					// three outputs
   // G = [1 0   (D/(1+D^3)
   //      0 1   (D^2/(1+D^3)]
   unsigned int g3num[] = {4,2};	// third column of G: D and D^2
   unsigned int g3den = 9;		    // =1001 denominator  1+D^3
   int p3 = 3;				// degree of denominators
   BinConvIIR conv3(k3,n3,p3,g3num,g3den);
   int in3 = 1;
   cout << "Third Test: " << endl;
   printconvtest(conv3,in3, Nimpulse);
   cout << endl << endl << endl;

}

void printconvtest(BinConv &conv, int in1, int Nimpulse)
{
   int i,j,j1;

   // test the getstate/setstate functions
   for(i = 0; i < (1<<conv.nu); i ++) {
	  conv.setstate(i);
	  cout << "i=" << i << "  getstate=" << conv.getstate() << endl;
   }
   cout << endl;

   // set up the impulse response inputs
   unsigned char *impulsein = new unsigned char[conv.k];
   unsigned int *out;
   cout << "Impulse reponse" << endl;
   conv.setstate(0);			// set back to state 0
   for(i = 0; i < Nimpulse; i++) {
	  if(i == 0) {				// set a pulse in indicated time slots
		 for(j1 = 0; j1 < conv.k; j1++) impulsein[j1] = (in1& (1<<j1)) != 0;
	  }
	  out = conv.encode(impulsein);
	  cout << "in: ";
	  for(j = 0; j < conv.k; j++) {
		 cout << int(impulsein[j]) << " ";
	  }
	  cout << "Out: ";
	  for(j = 0; j < conv.n; j++) {
		 cout << int(out[j]) << " ";
	  }
	  cout << "  State: " << conv.getstate() << endl;
	  for(j1 = 0; j1 < conv.k; j1++) impulsein[j1] = 0;
   }

   // print the state/nextstate table
   cout << "State/Next State Table" << endl;
   unsigned char in[conv.k];
   for(i = 0; i < (1<<conv.nu); i++) {  // for each state
	  cout << "state " << i << ": ";
	  for(j = 0; j < (1<<conv.k); j++) {	// for each input
		 conv.setstate(i);		// set the state
		 // convert the input as a number to a binary vector
		 for(j1 = 0; j1 < conv.k; j1++) in[j1] = (j & (1<<j1))!= 0;
		 out = conv.encode(in);		// apply the input
		 cout << conv.getstate() << " (";
		 for(j1 = 0; j1 < conv.k; j1++) cout << int(in[j1]);
		 cout << "/";
		 for(j1 = 0; j1 < conv.n; j1++) {  // also print the output
			cout << int(out[j1]);
		 }
		 cout << ") ";
	  }
	  cout << endl;
   }
}

/*
Local Variables:
compile-command: "g++ -o testconvenc -g testconvenc.cc BinConvFIR.cc BinConvIIR.cc"
End:
*/



//  Program: finddfree.cc --  find d_free for a given set of 
//   connection coefficients.
//  (Currently works only for k=1)
//  Todd K. Moon

#include <iostream>
using namespace std;

#include "BinConvFIR.h"
#include "BinConvdec01.h"

#include "matalloc.h"

long int octconv(char *str, int & deg1);
void printconvtest(BinConv &conv, int in1, int Nimpulse);
void leftadjfix(char *str);

int main()
{
   int n,k,n1,k1;
   int dfree;

   cout << "enter k: ";
   cin >> k;
   if(k != 1) {
	  cout << "Warning: currently only implemented for k=1\n";
   }
   cout << "enter n: ";
   cin >> n;
   unsigned int **g;
   CALLOCMATRIX(g,unsigned int, k,n);
   int *p1 = new int[k];						// degrees of rows
   char str[20];
   int deg1,deg;
   int gint;
   
   k = 1;
   k1 = 0;
   deg = 0;
   for(n1 = 0; n1 < n; n1++) {
	  cout << "Enter g (as an octal number): ";
	  cin >> str;
	  if(str[0] == '0' || str[0] == '1' || str[0] == '2' || str[0] == '3') {
		 leftadjfix(str);
		 cout << str <<  endl;
	  }
	  gint = octconv(str,deg1);
	  cout << "gint=" << gint << "  deg1=" << deg1 << endl;
	  g[k1][n1] = gint;
	  if(deg1 > deg) {
		 deg = deg1;
	  }
   }
   p1[k1] = deg;
   BinConvFIR conv1(k,n,p1,g);

   // printconvtest(conv1, 1, 7);

   cout << "nu=" << conv1.nu << endl;
   BinConvdec01 decoder(conv1,8);
   dfree = decoder.finddfree();
   cout << "dfree=" << dfree << endl;

}

long int octconv(char *str, int & deg1)
{
   int i;
   long int num = 0;
   long int num1 = 0;
   long int mask, mask1;

   deg1 = 0;
   for(i = 0; i < strlen(str); i++) {
	  num = num*8 + (str[i]-'0');
   }
   // determine the degree of the connection polynomial
   // assumes that the polynomial is "left-justified"
   mask = (1 << (strlen(str)*3-2));
   i = 1;
   while(mask){
	  if(mask & num) {
		 deg1 = i;
	  }
	  mask /= 2;
	  i++;
   }
   // flip the bits around
   num1 = 1;
   mask = (1 << (strlen(str)*3-2));
   mask1 = 2;
   while(mask) {
	  if(mask & num)
		 num1 |= mask1;
	  mask /= 2;
	  mask1 *= 2;
   }
   
   return num1;
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
   unsigned int nextst;

   // set up the impulse response inputs
   unsigned char *impulsein = new unsigned char[conv.k];
   unsigned char *out;
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
	  if(i==128) {
		 cout << "i=128\n";
	  }
	  for(j = 0; j < (1<<conv.k); j++) {	// for each input
		 conv.setstate(i);		// set the state
		 // convert the input as a number to a binary vector
		 for(j1 = 0; j1 < conv.k; j1++) in[j1] = (j & (1<<j1))!= 0;
		 out = conv.encode(in);		// apply the input
		 nextst = conv.getstate();
		 cout << nextst << " (";
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

void leftadjfix(char *str)
{
   int num = 0;
   int nshift = 0;
   int i, j,num2;
   char revstr[200];

   for(i = 0; i < strlen(str); i++) {
	  if(i == 0) {
		 num2 = num = str[i]-'0';
		 while(num2 < 4) {
			num2 <<= 1;
			nshift++;
		 }
	  }
	  else {
		 num = num*8 + (str[i] -'0');
	  }
   }
   num <<= nshift;			/* shift over */
   i = 0;
   while(num) {
	  revstr[i++] = (num % 8) + '0';
	  num /= 8;
   }
   revstr[i] = 0;
   str[i] = 0;
   for(j = i-1,i=0; j>= 0; j--,i++) {
	  str[j] = revstr[i];
   }
   printf("g0=%s\n",str);
}


/*
Local Variables:
compile-command: "g++ -o finddfree -Wno-deprecated -g finddfree.cc Convdec.cc BinConvFIR.cc BinConvdec01.cc"
End:
*/



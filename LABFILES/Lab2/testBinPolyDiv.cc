//  Program: BinLFSR.cc
//
//  Todd K. Moon
//

// #include <ostream.h>

extern "C" {
#include <stdio.h>
}

#include "BinPolyDiv.h"

int main()
{
   unsigned int g = 0x23;				// divisor: 10 0011 = 1+D+D^5
   unsigned char d[] = {1,1,0,0,0,1,0,1,1};	// dividend: 1+D+D^5 + D^7 + D^8
   unsigned char q[4];		  // room to hold a polynomial of degree 3
   int ddegree=8;
   int qdegree, rdegree;
   unsigned int rem;
   int remdeg;
   BinPolyDiv X(g,5);

   rem = X.div(d, ddegree, q, qdegree, rdegree);
   printf("rem=%x   rdegree=%d  qdegree=%d\n",rem,rdegree,qdegree);
   for(int i = 0; i <= qdegree; i++) {
	  printf("%d ",q[i]);
   }
   printf("\n");
 
}


/*
Local Variables:
compile-command: "g++ -o testBinPolyDiv -g testBinPolyDiv.cc BinPolyDiv.cc"
End:
*/



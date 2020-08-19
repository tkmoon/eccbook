//
//
//  Program: tcmrot1.cc --- test the constellation for 
//  a rotationally invariant code
//
//
//  Todd K. Moon
//  Utah State University
//
//  Date:  March 10, 2004

extern "C" {
#include <stdio.h>
}

int main()
{
   int nu = 3;					// number of bits of storage
   int st[nu];					// 
   int Nstate = (1<<nu);		// number of states
   int k = 2;					// number of input bits
   int Ninput = (1<<k);			// number of input values
   int ina[k];
   int out;
   int a,b,c,d;
   int state,input, mask, i,temp;

   for(state = 0; state < Nstate; state++) {
	  // write the state into the state array
	  for(mask = 1,i = 0; i < nu; i++,mask *= 2) {
		 if(mask & state) st[i] = 1; else st[i] = 0;
	  }
	  printf("state=(%d%d%d)\n",st[2],st[1],st[0]);
	  for(input = 0; input < Ninput; input++) {
		 // write the state into the state array
		 for(mask = 1,i = 0; i < nu; i++,mask *= 2) {
			if(mask & state) st[i] = 1; else st[i] = 0;
		 }
		 // write the input into the input array
		 for(mask=1,i=0; i < k; i++, mask *= 2) {
			if(mask & input) ina[i] = 1; else ina[i] = 0;
		 }
		 // compute at various points in the schematic
		 b = (st[1] + ina[1]) & 1;
		 c = (b & st[0]);
		 a = (ina[0] + ina[1] + c + st[2]) & 1;
		 d = ina[0] & st[0];
		 out = st[0];
		 temp = st[0];
		 // compute next state
		 st[0] = (b + d) & 1;
		 st[1] = a;
		 st[2] = temp;
		 printf("   in=(%d%d)  a=%d  b=%d  out=(%d%d%d)  nextstate=(%d%d%d)\n",
				ina[1],ina[0],a,b,ina[1],ina[0],out,st[2],st[1],st[0]);
	  }
   }
}
				
		 


/*
Local Variables:
compile-command: "g++ -o  tcmrot2 -Wno-deprecated -g tcmrot2.cc"
End:
*/



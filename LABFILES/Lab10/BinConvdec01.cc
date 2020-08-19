//  BinConvdec01.cc 
//  Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "BinConvdec01.h"
#include "matalloc.h"

void 
BinConvdec01::buildoutputmat(BinConv & encoder)
// builds the lookup from [state][input] to output array,
{
   unsigned int savestate = encoder.getstate();

   // outputmat[state][input][outputnum]
   CALLOCTENSOR(outputmat,unsigned char,numstates,numbranches,n);

   unsigned char ins[k];
   unsigned int state;
   unsigned int inp;
   int i,j;
   unsigned int *out;
   unsigned int outint;

   for(state = 0; state < numstates; state++) { // for each state
	  for(inp = 0; inp < numbranches; inp++)  {
		 encoder.setstate(state);
		 // convert inp to array
		 for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
		 out = encoder.encode(ins);
		 for(j = 0; j < n; j++) { 
			outputmat[state][inp][j] = out[j];
		 }
	  }
   }
   // cout << "outputs: " << endl;
   for(state = 0; state < numstates; state++) {
	  // cout << "state=" << state << ": ";
	  for(inp=0; inp < numbranches; inp++) {
		 for(j = 0; j < n; j++) {
			// cout << int(outputmat[state][inp][j]);
		 }
		 // cout << " ";
	  }
	  // cout << endl;
   }
   encoder.setstate(savestate);
}		 


/*
Local Variables:
compile-command: "g++ -c BinConvdec01.cc"
End:
*/

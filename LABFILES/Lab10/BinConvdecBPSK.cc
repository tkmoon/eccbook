//  BinConvdec01.cc 
//  Todd K. Moon

#include "BinConvdecBPSK.h"
#include "matalloc.h"
#include "BPSKmodvec.h"

void 
BinConvdecBPSK::buildoutputmat(BinConv & encoder)
// builds the lookup from [state][input] to output array,
{
   unsigned int savestate = encoder.getstate();

   // outputmat[state][input][outputnum]
   CALLOCTENSOR(outputmat,double,numstates,numbranches,n);

   unsigned char ins[k];
   unsigned int state;
   unsigned int inp;
   int i,j;
   unsigned int *out;
   unsigned int outint;
   double *modout;
   BPSKmodvec modulator(n);

   for(state = 0; state < numstates; state++) { // for each state
	  for(inp = 0; inp < numbranches; inp++)  {
		 encoder.setstate(state);
		 // convert inp to array
		 for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
		 out = encoder.encode(ins);
		 modout = modulator.mod(out);
		 for(j = 0; j < n; j++) { 
			outputmat[state][inp][j] = modout[j];
		 }
	  }
   }
   cout << "In BinconvdecBPSK: outputs: " << endl;
   for(state = 0; state < numstates; state++) {
	  cout << "state=" << state << ": ";
	  for(inp=0; inp < numbranches; inp++) {
		 for(j = 0; j < n; j++) {
			cout << int(outputmat[state][inp][j]);
		 }
		 cout << " ";
	  }
	  cout << endl;
   }
   encoder.setstate(savestate);
}		 


/*
Local Variables:
compile-command: "g++ -c BinConvdecBPSK.cc"
End:
*/

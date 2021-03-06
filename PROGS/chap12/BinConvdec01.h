// BinConvdec01.h -- Convolutional decoder for binary (0,1) data
// Todd K. Moon

#ifndef BinConvdec01_H
#define BinConvdec01_H
#include "Convdec.h"
#include "matalloc.h"

class BinConvdec01 : public Convdec {
   unsigned char ***outputmat;	// [numstates][numbranches][n] -- outputs
   unsigned char *data;			// pointer to input data
public:
   BinConvdec01(BinConv & encoder, int in_pathmem) 
	  : Convdec(encoder,in_pathmem)
	  {  
		 buildoutputmat(encoder);		// set up the outputmat
	  };
   ~BinConvdec01() { FREETENSOR(outputmat,numstates);  delete[] data;   }

   virtual double metric(unsigned int state, int branch) {
	  double sum=0;
	  for(int i = 0; i < n; i++) {
		 sum += (data[i] != outputmat[state][branch][i]);
	  }
	  return sum;
   };

   virtual double metric_dfree(unsigned int state, int branch)
   {
	  double sum=0;
	  for(int i = 0; i < n; i++) {
		 sum += outputmat[state][branch][i];
	  }
	  return sum;
   };


   int decode(unsigned char *outs) {
	  data = outs;
	  return viterbi();
   };
   int finddfree() {
	  return viterbifinddfree();
   };
   void buildoutputmat(BinConv& encoder);

};

#endif

/*
Local Variables:
compile-command: "g++ -c BinConvdec01.cc"
End:
*/

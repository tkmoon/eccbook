// BinConvFir.h -- declarations for an (n,k) polynomial
// binary convolutional coder
// Todd K. Moon
#ifndef BINCONVFIR_H
#define BINCONVFIR_H
#include "BinConv.h"
#include "matalloc.h"			// matrix allocation stuff

// The code is represented as a transfer function matrix
// G = [h11 h12 ... h1n
//      h21 h22 ... h2n
//      ...
//      hk1 hk2 ... hkn]
// 
// Each column of G must be an FIR filter As a result, each column is
// specified by a set of h vectors (the numerators)

#include "matalloc.h"

class BinConvFIR : public BinConv {
   unsigned int **h;			// encoder transfer functions (h0=lsb)
   int *nui;					// max row degree
   unsigned int *mem;			// state memory of encoders
   int maxdeg;					// maximum degree of a row (used by encode)
public:
   BinConvFIR(int in_k, int in_n, int *degs, unsigned int** h_in);
   ~BinConvFIR() { FREEMATRIX(h); delete[] nui; delete[] mem;
   };
   virtual unsigned int *encode(const unsigned char *ins); 
                                // encode one step of input
   virtual unsigned int getstate() const;
                                // return the state of the encoder
   virtual void setstate(const unsigned int state); 
                                // set the state of the encoder
};

#endif

/*
Local Variables:
compile-command: "gcc -c -g BinConvFIR.cc"
End:
*/

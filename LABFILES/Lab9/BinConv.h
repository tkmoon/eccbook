// BinConv.h -- declarations for the base class of a Binary 
// Convolutional encoder
// This is then specialized by inheritance to FIR or IIR encoders

// Todd K. Moon
#ifndef BINCONV_H
#define BINCONV_H

class BinConv {
protected:
   unsigned int *outs;			// array of n outputs
public:
   int n;						// number of outputs
   int k;						// number of inputs
   int nu;						// total of row degrees (constraint length)
   BinConv(int in_k, int in_n) {
	  k = in_k;
	  n = in_n;
	  outs = new unsigned int[n]; // output array
   }
   ~BinConv() { delete[] outs; };
   virtual unsigned int *encode(const unsigned char *ins) = 0;
                               // encode one step of input (pure virtual)
   virtual unsigned int getstate() const = 0;
                               // Get the state of the encoder
   virtual void setstate(const unsigned int state) = 0;
                               // Set the state of the encoder
};

#endif

/*
Local Variables:
compile-command: "gcc -c -g BinConv.cc"
End:
*/

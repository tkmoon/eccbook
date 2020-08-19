// BinTF2.h -- declarations for Multi-input Binary Transfer function class
// Todd K. Moon
#ifndef BINTF2_H
#define BINTF2_H

// Implement the multi-input relationship
// y(z) = H1(z) in1(z) + H2(z) in2(z) + ... Hk(z) ink(z)
// where all Hi share a common denominator,
// Hi(z) = Pi(z)/Q(z)
// This can be represented in matrix form:
// y(z) = [in1(z) .. ink(z)] [  P1(z)/Q(z)
//                               ...
//                              Pk(z)/Q(z)
// 
//     h0 z^{-p} + h1 z^{-p+1} + ... + hp      h0 + h1 z + ... + hp z^p
//H(z)=----------------------------------  =   ------------------------
//     g0 z^{-p} + g1 z^{-p+1} + ... + gp       g0 + g1 z + ... + gp z^p
//
// over GF2, using bits of an integer to represent the state
// (so p < 32)
// The feedforward connection is represented using two integers,
// h and g, whose  LSB are h0 and g0, resp.
//
// For example: H(z) = 1+z^{-1}: p=1,   h0=1,h1=1  -> h=3
//                                      g0=0,g1=1  -> g=2
// However, in the FIR case it suffices to set g=0 (to indicate FIR filter)
//
// For example, H(z) = z/(z+1) : p=1,   h0=0 h1=1  -> h=2
//                                      g0=1 g1=1  -> g=3
//
// For an FIR filter, use g=0 (default argument)
//
// This can be somewhat confusing for IIR, because things are often expressed
// in terms of D, not z:
// Consider the transfer function  
//     D/(1+D^3) = z^(-1)/(1+z^(-3)) = z^2/(z^3+1)
// In this case, h0=h1=0  h2=2  and g0=1, g1=0  g2=1
// so h = 4  and g = 5.  One might have thought (without some care)
// that h=2.
//
// For D^2/(1+D^3) = z/(z^3+1) we have h=2 and g = 5.
// So, Be Careful!


class BinTF2 {
   int p;						// denominator degree
   int k;						// number of transfer functions in system
   unsigned int state;			// state of system
   unsigned int mask;
   unsigned int g;				// denominator
   unsigned int *h;				// numerators
   unsigned int *hshift;			// shifted numerator
   int pm1;						// p-1;
public:
   BinTF2(void) { // default constructor
	  p=k=state=mask=g=0;
	  h=hshift=0;
   };
   // real constructor
   BinTF2(int in_p, int in_k, unsigned int *h_in, unsigned int g_in=0) {
	  p = in_p;
	  pm1 = p? p-1 : 0;
	  k = in_k;
	  h = new unsigned int[k];
	  hshift = new unsigned int[k];
	  mask = (1<<p)-1;
	  g = g_in & mask;
	  for(int i = 0; i < k; i++) {
		 h[i] = h_in[i];
		 hshift[i] = p? (h[i]? (h[i]>>1):0) : h[i];
	  }
	  state = 0;
   }
   // "build" the object according the to passed-in parameters
   void set(int in_p, int in_k, const unsigned int *h_in, unsigned int g_in=0)
   {
	  if(h) delete[] h;
	  if(hshift) delete[] hshift;
	  p = in_p;
	  pm1 = p? p-1 : 0;
	  k = in_k;
	  h = new unsigned int[k];
	  hshift = new unsigned int[k];
	  mask = (1<<p)-1;
	  g = g_in & mask;
	  for(int i = 0; i < k; i++) {
		 h[i] = h_in[i];
		 hshift[i] = p ? (h[i]? (h[i]>>1) : 0) : h[i];
	  }
	  state = 0;
   };
   // Given an input array of length k, step the transfer function object
   unsigned int step(const unsigned char *in) {
	  int i;
	  unsigned int out = state;
	  for(i = 0; i < k; i++) {
		 out ^=  hshift[i]*in[i];
	  }
	  out >>= pm1;
	  unsigned int st = (state<<1)^(g*out);
	  for(i = 0; i < k; i++) {
		 st ^= h[i]*in[i];
	  }
	  state = st&mask;
	  return out;
   }
   // get the transfer function state
   unsigned int getstate() const { return state; }
   // set the transfer funtion state
   void setstate(unsigned int instate) { state = instate; }
   // interface for internal information
   unsigned int getp() const { return p;}
   unsigned int getg() const { return g;}
};

#endif

/*
Local Variables:
compile-command: ""
End:
*/

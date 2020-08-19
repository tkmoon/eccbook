// Convdec.h --- Convolutional decoding using the Viterbi algorithm
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef Convdec_H
#define Convdec_H

#include "BinConv.h"
#include <iostream>
using namespace std;

class Convdec {
protected:
   int n;						// number of outputs
   int k;						// number of inputs
   int numstates;				// number of states
   int numbranches;				// number of branches (2^k)
   int viterbi();				// decode using the data int data
   int viterbifinddfree();		// find the free distance using the VA 


   unsigned int **prevstate;	// [numstates][numbranches] -- states prior
   unsigned int **inputfrom;	// [numstates][numbranches] -- 
                                // inputs from prior states

   double *metrics1;			// path metrics to each state
   double *metrics2;			// path metrics to each state
   // (metrics1 and metrics2 are used to double buffer and avoid copying)
   double *metrics;				// point to either metric1 or metric2
   double *othermetrics;		// used to swap

   unsigned int startstate;		// the starting state of trellis (default=0)
   void buildprev(BinConv &encoder); // built array of previous states
   virtual void buildoutputmat(BinConv & encoder)=0; // built output matrix

   // Viterbi Variables:
   int pathmem; // path length

   struct VITPATH {
	  int *state;
	  int *input;
   };

   VITPATH *paths; // the list of paths leading to the ith time instant

   // paths[i].state[state] tells what state LEADS to state at the ith branch
   // paths[i].input[state] tells what input got to state from the vector
   // previous state
   
   // Thus, in a sense, paths[].state leads backward, and paths[].input
   // leads forward.

   // paths is stored in as a circular queue, with pointers fpath and bpath

   // used int the algorithm:
   int fpath; // the pointer to the front of the paths
   int bpath; // the pointer to the back of the path
   int numbranchdec;   // number of branches decoded  -- count to keep
								// track of whether this is greater than the
								// desired memory


   void setpaths();				// set default path values
   
   virtual double metric(unsigned int state, int branch)=0;
   // compute path metric
   virtual double metric_dfree(unsigned int state, int branch) = 0;

   void setstartstate(unsigned int ststate) {
	  startstate=ststate;
	  setpaths();
   }

   int inci(int i) { // increment a counter through the circular path list
	  i = (i+1) % pathmem;
	  return(i);
   }

   int deci(int i) { // decrement a counter through the circular path list
	  i = i-1;
	  if(i < 0) i+= pathmem;
	  return(i);
   }
   unsigned int beststate;		// best current state
public:

   Convdec(BinConv &encoder, int in_pathmem);
   // return 1 if decoded value is available; inputs are written into 'inputs'
   void showpaths(void);		// dump the paths out for debugging
   unsigned int inputs;	  // the inputs found by the Viterbi algorithm
                          // (as a bit array)
   int getinpnow(int adv);	// Get inputs from current trellis
                            // Return 1 if there is new input available
                            // If adv is set, then advance pointer in path
                            // Write input value into 'inputs'

};


#endif
/*
Local Variables:
compile-command: "g++ -c -g Convdec.cc"
End:
*/

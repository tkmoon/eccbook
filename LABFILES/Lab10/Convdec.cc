//  Program: Convdec.cc -- Convolutional decoding (using the Viterbi alg)
//  Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "Convdec.h"
#include "BinConv.h"
#include "matalloc.h"
#include <iostream>
#define LARGE 1e99
#include <iostream>

Convdec::Convdec(BinConv &encoder, int in_pathmem)
{
   n = encoder.n;
   k = encoder.k;

   pathmem = in_pathmem;

   // build up the information needed to do the decoding
   numstates = (1<<encoder.nu);
   numbranches = (1<<k);
   metrics1 = new double[numstates];
   metrics2 = new double[numstates];

   buildprev(encoder);   // build the state/previous state array

   // build up the path information
   paths = new VITPATH[pathmem];
   for(int i = 0; i < pathmem; i++) {
	  paths[i].state = new int[numstates];
	  paths[i].input = new int[numstates];
   }
   startstate = 0;				// default state state
   setpaths();					// initialize the path variables
}

int
Convdec::viterbi()
// r is an array of the received values
{
   int i;
   double d = LARGE;			// metric to best state
   unsigned int state;			// state index variable
   double bestcost;				// best cost to a state
   double cost;					// path cost to state
   int branch;					// branch number from state
   int bestinp;					// the best input from state
   unsigned int fromst;			// from state
   unsigned int bestfrom;		// the best from state
   unsigned int bestst_back;	// the best state as we work back
   double *metswap;				// used to swap metric arrays
   unsigned int fromin;			// input from previous state


   // Fill in the blanks ...

   return(1);
}

int
Convdec::getinpnow(int adv)
// Get the best input from the current state of the trellis
// if adv is set, move forward
{
   // if we reach here, we are ready to decode
   // we already know the best one, from the assignment above.  
   // Only need to follow it back
   if(bpath == fpath) return 0;  // no information to print

   unsigned int bestst_back;
   int i = deci(fpath);
   bestst_back = beststate;
   do {
	  bestst_back = paths[i].state[bestst_back];
	  i = deci(i);
   }
   while(i != bpath);
   inputs = paths[bpath].input[bestst_back];
   if(adv) {
	  bpath = inci(bpath);
   }
   return 1;
}

void 
Convdec::buildprev(BinConv &encoder)
{
   unsigned int savestate = encoder.getstate();

   CALLOCMATRIX(prevstate,unsigned int, numstates, numbranches);
   CALLOCMATRIX(inputfrom,unsigned int, numstates, numbranches);

   // first build the nextstate 
   unsigned int **nextstate;
   CALLOCMATRIX(nextstate,unsigned int,numstates,numbranches);
   unsigned char ins[k];
   unsigned int state;
   unsigned int inp;
   int i;
   unsigned int nextst;
   int *nfrom = new int[numstates];

   for(state = 0; state < numstates; state++) {
	  for(inp = 0; inp < numbranches; inp++)  {
		 encoder.setstate(state);
		 // convert inp to array
		 for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
		 encoder.encode(ins);
		 nextst = encoder.getstate();
		 prevstate[nextst][nfrom[nextst]] = state;
		 inputfrom[nextst][nfrom[nextst]] = inp;
		 nfrom[nextst]++;
	  }
   }
	  
   for(state = 0; state < numstates; state++) { // for each state
	  for(inp = 0; inp < numbranches; inp++)  {
		 encoder.setstate(state);
		 // convert inp to array
		 for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
		 encoder.encode(ins);
		 nextstate[state][inp] = encoder.getstate();
	  }
   }
   for(state = 0; state < numstates; state++) {
	  unsigned int *outs;
	  // cout << "state=" << state << ": ";
	  for(inp=0; inp < numbranches; inp++) {
		 // cout << nextstate[state][inp] << " ";
		 encoder.setstate(state);
		 for(i = 0; i < k; i++) { if(inp&(1<<i)) ins[i] = 1; else ins[i] = 0;}
		 outs = encoder.encode(ins);
		 // cout << "(";
		 for(i = 0; i < n; i++) {
			// cout << int(outs[i]);
		 }
		 // cout << ") ";
	  }
	  // cout << endl;
   }
   encoder.setstate(savestate);


   // print the prevstate table
   // cout << "fromstates: " << endl;
   for(state = 0; state < numstates; state++) {
	  // cout << "state=" << state << ": ";
	  for(inp=0; inp < (1<<k); inp++) {
		 // cout << prevstate[state][inp] << " (";
		 // cout << inputfrom[state][inp] << ") ";
	  }
	  // cout << endl;
   }
   delete[] nfrom;
   encoder.setstate(savestate);
   FREEMATRIX(nextstate);
}



void Convdec::setpaths()
{
   int i;

   fpath = 0;
   bpath = 0;
   numbranchdec = 1;
   for(i = 0; i < numstates; i++) {
	  metrics1[i] = LARGE;
   }
   metrics1[startstate] = 0;
   metrics = metrics1;
   othermetrics = metrics2;
}


void Convdec::showpaths(void)
{
   int i,state;

   return;  // don't print ....

   cout << "fpath=" << fpath << "\tbpath=" << bpath << "\n";

   cout << "Metrics: \n";
   for(state = 0; state < numstates; state++) {
	  cout << metrics[state] << "\t";
   }


   cout << "\n\nPrevious State Path\n";
   for(state = 0; state < numstates; state++) {
	  for(i = bpath; i != fpath; i = inci(i)) {
		 cout <<paths[i].state[state] << "\t";
	  }
	  cout << "\n";
   }

   cout << "\nInput path\n";
   for(state = 0; state < numstates; state++) {
	  	  for(i = bpath; i != fpath; i = inci(i)) {
		 cout <<paths[i].input[state] << "\t";
	  }
	  cout << "\n";
   }
}


int
Convdec::viterbifinddfree()
// r is an array of the received values
{
   int i;
   double d = LARGE;			// metric to best state
   unsigned int state;			// state index variable
   double bestcost;				// best cost to a state
   double smallestcost;			// smallest cost to any state
   double cost;					// path cost to state
   int branch;					// branch number from state
   int bestinp;					// the best input from state
   unsigned int fromst;			// from state
   unsigned int bestfrom;		// the best from state
   unsigned int bestst_back;	// the best state as we work back
   double *metswap;				// used to swap metric arrays
   unsigned int fromin;			// input from previous state

   int nbranch = 0;
   int dfree = 65535;
   int first = 1;
   do {
	  nbranch++;
	  smallestcost = LARGE;
//cout << "nbranch=" << nbranch << endl;
//cout << "path costs: ";
	  for(state = 0; state < numstates; state++) { 
		 // assume that we have data fromstates[state][in], which provides
		 // a list of all 'from' states to 'state'
		 bestcost = LARGE;
		 for(branch = 0; branch < numbranches; branch++) {
			fromst = prevstate[state][branch];
			if(state == 0 && fromst == 0) continue; // 
			if(fromst == 0 && nbranch > 1) continue;
			fromin = inputfrom[state][branch];
			cost = metrics[fromst] + metric_dfree(fromst, fromin);
// cout << "state=" << state << " fromst=" << fromst << " cost=" <<
//			   metric_dfree(fromst,fromin) << " ";
			if(cost < bestcost) {
			   bestcost = cost;
			   bestfrom = fromst;
			   bestinp = fromin;
			}
		 }
		 // record path information
		 othermetrics[state] = bestcost;
		 if(bestcost < smallestcost) { //  find the smallest cost to all states
			smallestcost = bestcost;
		 }
//cout << othermetrics[state] << " ";
	  }
	  // swap the metric information
	  metswap = metrics;   metrics = othermetrics;  othermetrics = metswap;

	  if(metrics[0] < dfree) {		// we are back to state 0
		 dfree = int(metrics[0]);
	  }
//cout << "   dfree=" << dfree << "  smallestcost=" << smallestcost;
//cout << endl;
	  if(smallestcost >= dfree) { // no other branches will be smaller
		 break;
	  }
   } while(1);
   return dfree;
}

/*
Local Variables:
compile-command: "g++ -c -g Convdec.cc"
End:
*/



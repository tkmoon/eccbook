// BinLFSR.h -- declarations for BinLFSR class
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef BinLFSR_H
#define BinLFSR_H

class BinLFSR {
   unsigned int g;
   int state;
   unsigned int mask;
   unsigned int mask1;
   int n;
public:
   BinLFSR(void) { g=n=state=mask=mask1=0;}// default constructor
   BinLFSR(unsigned int g, int n, unsigned int initstate=1);// constructor
   ~BinLFSR() {};                        // destructor
   void setstate(int state);             // set the state
   char step(void);                      // step once; return output
   char step(unsigned int &state);       // step once; return output and state
   void steps(int nstep, char *outputs); // step nstep times
};
#endif

/*
Local Variables:
compile-command: "gcc -c -g BinLFSR.cc"
End:
*/

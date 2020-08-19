//  Program: BinLFSR.cc
//
//  Todd K. Moon
//

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "BinLFSR.h"

BinLFSR::BinLFSR(unsigned int in_g, int in_n, unsigned int initstate)
{
   int i;
   n = in_n;
   mask = (1<<n)-1;  // a mask of n ones
   mask1 = (1<<(n-1));  // select the highest bit in shift register
   g = in_g & mask;
   state = initstate;
}   

char BinLFSR::step(void)
{
   char out = (state&mask) >> (n-1);
   state = ((state<<1)^(g*out))&mask;
   return out;
}

char BinLFSR::step(unsigned int &stateout)
{
   char out = (state&mask) >> (n-1);
   stateout = state;
   state = ((state<<1)^(g*out))&mask;
   return out;
}
/*
Local Variables:
compile-command: "gcc -c -g BinLFSR.cc"
End:
*/



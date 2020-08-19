//  Program: BinLFSR.cc
//
//  Todd K. Moon
//
// Copyright 2019 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include <iostream>
using namespace std;

extern "C" {
#include <stdio.h>
}


#include "BinLFSR.h"

int main()
{
   unsigned int g = 0x17; // 10111 = 1+d+d^2+d^4
   // int g = 0x13; // 10011 = 1+d+d^4
   int Nsteps = 16;
   BinLFSR X(g,4);
   int i;
   char out;
   unsigned int state;

   for(i = 0; i < Nsteps; i++) {
	  out = X.step(state);
	  cout << out << " " << state << "\n";
   }
}


/*
Local Variables:
compile-command: "g++ -o testBinLFSR -g testBinLFSR.cc BinLFSR.cc"
End:
*/



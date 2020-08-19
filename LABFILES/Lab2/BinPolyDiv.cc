// BinPolyDiv.h -- declarations for BinPolyDiv class
// Todd K. Moon

#include "BinPolyDiv.h"

BinPolyDiv::BinPolyDiv(unsigned int in_g, int in_p)
{
   p = in_p;
   mask = (1<<p)-1;
   g = in_g & mask;
   state = 0;
}

unsigned int BinPolyDiv::div(unsigned char *d, int n, unsigned char *q, 
					int &quotientdegree, int &remainderdegree)
// divide d by g. Return quotient in q and remainder as return value (integer)
// d has d0 first, so start at upper end of d and work back.
{
   int i,j;
   unsigned char out;
   int outdeg=n-p;
   quotientdegree=0;
   state = 0;

   // Fill in the blanks ...
   return state;
}

unsigned int BinPolyDiv::remainder(unsigned char *d, int n,int &remainderdegree)
// divide d by g. Return remainder as return value (integer)
// d has d0 first, so start at upper end of d and work back.
{
   int i,j;
   unsigned char out;
   state = 0;

   // Fill in the blanks ...
   return state;
}

unsigned int BinPolyDiv::rem1(unsigned char *d, int n)
// divide d by g. Return remainder as return value (integer)
// d has d0 first, so start at upper end of d and work back.
// This function does not bother to determine the remainderdegree
{
   int i,j;
   unsigned char out;

   // Fill in the blanks ...
   return state;
}

unsigned int BinPolyDiv::rem2(unsigned char *d, int n,int &remainderdegree)
// divide d by g. Return remainder as return value (integer)
// d has d0 first, so start at upper end of d and work back.
// This function does not clear the state to 0, and does not perform 
// the initial shift,  so incremental operations can be performed.
{
   int i,j;
   unsigned char out;

   // Fill in the blanks ...
   return state;
}


/*
Local Variables:
compile-command: "gcc -c -g BinPolyDiv.cc"
End:
*/

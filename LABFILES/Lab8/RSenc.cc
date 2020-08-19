// RSenc.cc -- a general RS encoder
// Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "RSenc.h"
#include "polynomialT.cc"

template class polynomialT<GFNUM2m>;
template class polytemp<GFNUM2m>;

RSenc::RSenc(int n_in, int k_in, int t_in, int j0_in, GFNUM2m A1_in)
// constructor
{
   n = n_in;
   k = k_in;
   t = t_in;
   j0 = j0_in;
   A1 = A1_in;
   g = GFNUM2m(1);
   GFNUM2m lc[2] = {1,1};
   polynomialT<GFNUM2m> l(1,lc);  // set up the linear term
   for(int i=0; i < 2*t; i++) {
	  l[0] = (A1^(j0+i));
	  g *= l;
   }
   m.setc(k-1,1);				// allocate space for message polynomial
}

GFNUM2m *
RSenc::encode(GFNUM2m *m_in)
// encode from Galios field elements
{
   // fill in the blanks ...
}

void
RSenc::encode(unsigned char *m_in, unsigned char *c_out)
// encode from chars
{
   // fill in the blanks
}
   
/*
Local Variables:
compile-command: "g++ -c -g RSenc.cc"
End:
*/
	  

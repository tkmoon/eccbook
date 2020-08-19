// ChienSearch.h --- perform a Chien Search
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef ChienSearch_H
#define ChienSearch_H

#include "GFNUM2m.h"

class ChienSearch {
   int t;   // maximum degree
   GFNUM2m *Lambdas;  // registers for algorithm
   GFNUM2m *roots;    // stored root values
   GFNUM2m *alphapowers;  // powers of alpha (to speed up)
   int nroots;
public:
   ChienSearch(int in_t) {
	  t = in_t; nroots = 0;
	  roots = new GFNUM2m[t]; 
	  alphapowers = new GFNUM2m[t+1];
	  Lambdas = new GFNUM2m[t+1];
	  for(int i = 1; i <= t; i++) 
		 alphapowers[i] = A^i;
   }
   ~ChienSearch() { delete[] roots;  delete[] Lambdas; delete[] alphapowers;};
   GFNUM2m *getroots() { return roots; }
   int getnroots() { return nroots; }
   GFNUM2m *Search(GFNUM2m *poly,int nu, int & nroots);
                                // perform the Chien search
								//  on the polynomial of degree nu, 
                                //  returning a pointer to the roots array
                                // and the number of roots actually found
  
};

#endif
/*
Local Variables:
compile-command: "g++ -c -g ChienSearch.cc"
End:
*/

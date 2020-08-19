//
//
//  Program:  fht.cc
//  Compute the fast Hadamard transform on the integer sequence F
//  where F has 2^m points in it.
//
//  Todd K. Moon
//  Utah State University
//
//  Date:  March 24, 2004
//


void fht(int *F, int m)
{
   int i;
   int n;
   int skip,skip2;
   int j, j1, j1skip, k1;
   int tmp;

   n = (1<<m);					// number of points
   skip = 1;
   skip2 = 2;
   for(i = 0; i < m; i++) {		// over each iteration
	  j = 0;					// start at the first line
	  while(j < n) {
		 j1 = j;
		 j1skip = j1+skip;
		 for(k1 = 0, j1=j,j1skip=j1+skip; k1 < skip; k1++,j1++,j1skip++) {
			tmp = F[j1];
			F[j1] += F[j1skip];
			F[j1skip] = tmp - F[j1skip];
		 }
		 j += skip2;
	  }
	  skip = skip2;
	  skip2 <<= 1;
   }
}

void fhtinv(int *F, int m)
{
   int N = 1<<m;
   fht(F,m);
   for(int i = 0; i < (1<<m); i++) {
	  F[i] /= N;
   }
}


/*
Local Variables:
compile-command: "g++ -o testfht testfht.cc fht.cc"
End:
*/



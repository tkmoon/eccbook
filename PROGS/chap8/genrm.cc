//
//
//  Program: genrm.cc
//  Create a generator for an RM(r,m) code
// This is a kind of quick, brute-force (non recursive) way to 
// do this, that can only be coded explicitly up to a fixed value of r
// But it is easy to generalize to others, if desired.
//
//  Todd K. Moon
//  Utah State University
//
//  Date:  Jan 14, 2005
//  Last update:
//
#include <iostream>
using namespace std;
#include "matalloc.h"

int binom(int n, int k);

int main(int argc, char *argv[])
{
   int r, m;
   int n;
   int k;
   unsigned int i,i1,rc,j,i2,i3,i4;
   unsigned int shift,shift1,shift2,shift3,shift4;
   if(argc == 1) {
	  cout << "Usage: " << argv[0] << " r m" << endl;
	  exit(-1);
   }

   r = atoi(argv[1]);
   m = atoi(argv[2]);

   n = (1<<m);					// code length
   k = 1;
   for(i = 1; i <= r; i++) {		// compute the code dimension
	  k += binom(m,i);
   }
   cout << "k=" << k << endl;
   unsigned char **G;
   CALLOCMATRIX(G,unsigned char,k,n);
   for(int i = 0; i < k; i++) {
	  for(int j = 0; j < n; j++) {
		 G[i][j] = 0;
	  }
   }
   rc = 0;						// row count
   for(i = 0; i < n; i++) {		//  build the first row
	  G[rc][i] = 1;
   }
   if(r >= 1) {
	  for(i1 = 0; i1 < m; i1++) {	// order(1) operations
		 shift = (1<<(m-i1-1));
		 rc++;
		 for(i = 0; i < n; i++) {
			if(i & shift) {
			   G[rc][i] = 1;
			}
		 }
	  }
   }
   if(r >= 2) {
	  for(i1 = 0; i1 < m; i1++) {	// order(2) operations
		 shift1 = (1<<(m-i1-1));
		 for(i2 = i1+1; i2 < m; i2++) {
			if(i1==i2) continue;
			shift2 = (1<<(m-i2-1));
			rc++;
			for(i = 0; i < n; i++) {
			   if((i & shift1) && (i&shift2)) {
				  G[rc][i] = 1;
			   }
			}
		 }
	  }
   }

   if(r >= 3) {
	  for(i1 = 0; i1 < m; i1++) {	// order(3) operations
		 shift1 = (1<<(m-i1-1));
		 for(i2 = i1+1; i2 < m; i2++) {
			shift2 = (1<<(m-i2-1));
			for(i3 = i2+1; i3 < m; i3++) {
			   shift3 = (1<<(m-i3-1));
			   rc++;
			   for(i = 0; i < n; i++) {
				  if((i & shift1) && (i&shift2) && (i&shift3)) {
					 G[rc][i] = 1;
				  }
			   }
			}
		 }
	  }
   }

   if(r >= 4) {
	  for(i1 = 0; i1 < m; i1++) {	// order(4) operations
		 shift1 = (1<<(m-i1-1));
		 for(i2 = i1+1; i2 < m; i2++) {
			shift2 = (1<<(m-i2-1));
			for(i3 = i2+1; i3 < m; i3++) {
			   shift3 = (1<<(m-i3-1));
			   for(i4 = i3+1; i4 < m; i4++) {
				  shift4 = (1<<(m-i4-1));
				  rc++;
				  for(i = 0; i < n; i++) {
					 if((i & shift1) && (i&shift2) && (i&shift3)
						&& (i&shift4)) {
						G[rc][i] = 1;
					 }
				  }
			   }
			}
		 }
	  }
   }
   
   if(r >= 5) {
	  cout << "Not implmented for r in this range\n";
	  exit(0);
   }
   


   for(j = 0; j < k; j++) {
	  for(i = 0;i < n; i++) {
		 cout << int(G[j][i]);
	  }
	  cout << "\n";
   }
}
	  
   
   
int binom(int n, int k)
{
   double prod = 1;
   long int ipod;
   double i,j;
   if(k > n) return 0;
   if(k > n/2) k = n-k;
   if(k <= 1) {
	  if(k == 0) return 1;
	  if(k == 1) return n;
   }
   for(i = n-k+1, j=1; i <= n; i++,j++) {
	  prod *= i/j;
   }
   ipod = (long int)(prod+0.5);
   return ipod;
}


/*
Local Variables:
compile-command: "g++ -o genrm -g genrm.cc"
End:
*/



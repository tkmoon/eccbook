/**************************************
*
*  Program: primfind.c --
*  Program for generating binary primitive polynomials
*
*  This program was modified and commented from the program in
* Appendix A of Wicker 1993.
*
*  Todd K. Moon
*  Utah State University
*
*  Date: March 30, 2004
*
*  Example:
*  primfind -p 3 -s 3 -u 5 -w 3
*  Find primitive polynomials modulo 3 of weight 3 with degrees from
*   3 to 5
***************************************/

/* The program works by walking recursively over all d sequences (with the top 
   coefficient equal to 1 for convenience).  For each sequence, an 
   LFSR sequence is generated, starting in the all-1 state.  The LFSR
   is clocked until the initial state repeats, and the number of steps
   is counted.  If the number of steps is equal to p^n-1, then the 
   sequence represents a primitive polynomial.
*/

#include <stdlib.h>
#include <string.h>
#define N 32
#include <stdio.h>

/* the following data are globally used, and across
   recursive calls, so it makes sense for them to be declared global. */

int n;							/* degree of polynomial */
char d[N];						/* current set of coefficients */
char s[N];						/* state register for LFSR */
int p = 3;						/* modulus */

/* the following are local variables, but are declared here global
   so they don't have to be pushed on the stack for the recursive calls.
*/
char t;							/* LFSR out */
char f;							/* loop stop flag */
int i;							/* loop counter */
unsigned long c;				/* number of cycles counted before repeat */
unsigned long max;				/* max number of cycles before repeat */

void printd()  /* print out the array d */
{

  for(i = n; i >= 0; i--)
	printf("%c",d[i]+'0');
  printf("\n");
}

void LFSRgen()
{
   for(i = 0; i < n; i++) /* load the initial state of the LFSR */
	  s[i] = 1;
   c = 0;
   do {  /* repeat until initial state load repeats */
	  c++;
	  for(i = t = 0; i < n; i++) { /* propagate through a step of LFSR */
		 t = (t + s[i]*d[i]) % p;
	  }
	  for(i = 0; i < n-1; i++) { /* shift LFSR state */
		 s[i] = s[i+1];
	  }
	  s[n-1] = t;  /* feedback */
	  for(i = f = 0; i < n; i++) { 
		 if(s[i] != 1) { /*if any of these are not the initial load */
			f = 1;  /* keep going */
			break;
		 }
	  }
   }  while(f);
   /* if we reach here, then the LFSR has repeated.  If count==max,
      then we have a maximal length sequence */
   if(c == max) {
	 printd();  /* print out the current d */
   }
}

void gp(char l)
{
  int i;
  if(l==-1) {
	LFSRgen();	
  }
  else {
	if(l==0) {
	  for(i = 1; i < p; i++) {
		d[l] = i;
		gp(l-1);
	  }
	}
	else {
	  for(i = 0; i < p; i++) {
		d[l] = i;
		gp(l-1);
	  }
	}
  }
}

void gc(char l, char rw)
{
   char q;
   int i,j;
   if(rw==2) {
	  LFSRgen();
	  return;
   }
   for(q = l; q >= rw-2; q--) {
	 for(j = 1; j < p; j++) {	/* set the left digit */
	   d[0] = j;
	   for(i = 1; i < p; i++) {
		 d[q] = i;
		 gc(q-1,rw-1);
		 d[q] = 0;
	   }
	 }
   }
}

int main(int argc, char *argv[])
{
  int xx;
  int startn, endn;
  int w;

  if(argc == 1) {
	printf("\nUsage: %s [-p p] [-s s] [-u u] [-w w]\n\n",argv[0]);
	printf("Prints primitive polynomials modulo p\n");
	printf("-p p -- specifies the modulo.  Default: 2\n");
	printf("-s s -- specifies the starting degree.  Default: 2\n");
	printf("-u u -- specified the ending degree.  Default: 5\n");
	printf("-w w -- weight of polynomial.  Default: any\n");
	return(-1);
  }
  p = 2;
  startn = 2;
  endn = 5;
  w = 0;						/* any */
  for(i = 1; i < argc; i++) {
	if(!strcmp(argv[i],"-p"))
	  p = atoi(argv[++i]);
	if(!strcmp(argv[i],"-s"))
	  startn = atoi(argv[++i]);
	if(!strcmp(argv[i],"-u"))
	  endn = atoi(argv[++i]);
	if(!strcmp(argv[i],"-w"))
	  w = atoi(argv[++i]);
  }
  for(n = startn; n <= endn; n++) { /* this line determines the interval */
	printf("%d\n",n); 
	for(xx = max = 1; xx <= n; xx++) {
	  max = p*max;
	}
	max--;  /* max is p^n-1 */
	d[n] = 1;  /* set high-order coefficient */
	if(!w) {					/* select any weight */
	  gp(n-1);
	}
	else {
	  gc(n-1,w);
	}
	
	printf("\n");
  }  /* for n */
}



/*
Local Variables:
compile-command:"gcc -o primfind -g primfind.c"
End:
*/

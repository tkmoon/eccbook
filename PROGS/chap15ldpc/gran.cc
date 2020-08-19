/**********************************************************************/
/* Generate a unit-variance Gaussian random number.
   This function calls rand(), so you can control the seed using
   void srand(unsigned int seed);

   Following Numerical Recipes
*/
#include <stdlib.h>
#include <math.h>

static int graniset = 0;
static double grangset;

double gran(void)
{
   double rsq, v1, v2, fac;

   if(!graniset) {
	  graniset = 1;
	  do {
		 v1 = 2*(rand()/(double)RAND_MAX) - 1;
		 v2 = 2*(rand()/(double)RAND_MAX) - 1;
		 rsq = v1*v1 + v2*v2;
	  } while(rsq > 1 || rsq == 0);
	  fac = sqrt(-2*log(rsq)/rsq);
	  grangset = v1*fac;
	  return v2*fac;
   }
   else {
	  graniset = 0;
	  return grangset;
   }
}

void gran(double &r1, double &r2)
{
   double rsq, v1, v2, fac;

   do {
	  v1 = 2*(rand()/(double)RAND_MAX) - 1;
	  v2 = 2*(rand()/(double)RAND_MAX) - 1;
	  rsq = v1*v1 + v2*v2;
   } while(rsq > 1 || rsq == 0);
   fac = sqrt(-2*log(rsq)/rsq);
   r1 = v1*fac;
   r2 = v2*fac;
}

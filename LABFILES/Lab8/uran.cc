/* Generate a uniform random number between 0 and 1.  This 
   function calls rand(), so you can control the seed with srand().
*/

extern "C" {
#include <stdlib.h>
}

double uran(void)
{
   return rand()/(double)(RAND_MAX+1);
}

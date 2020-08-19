
// utiltkm.cc --
// some simple useful utilities 
// uran --- generates U(0,1) r.v.
// gran -- generators N(0,1) r.v.s, or pairs of them
// sort2lfd -- sort a double array, and a corresponding integer (index) array
// sort1lf -- sort a double array
// sort2fd -- sort a float array, and a corresponding integer (index) array
// sort1lf -- sort a float array

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


/* Generate a uniform random number between 0 and 1.  This 
   function calls rand(), so you can control the seed with srand().
*/


double uran(void)
{
   return rand()/(RAND_MAX+1.0);
}

int intran(int a, int b)
// Generate a uniform random number between a and b (inclusive)
{
   int r = int(a + floor(double(b-a)*rand()/(RAND_MAX+1.0)));
   return r;
}


/**********************************************************************
 *
 *  sort2.c -- a quicksort routine
 *
 **********************************************************************/

static double *alf;
static int *a2d;
static float *af;

/**********************************************************************/

static void quicksort2lfd(int left, int right)
{
    long i, j;
    double ref;
    int ref2;
    i = left;
    j = right;
    ref = alf[i];
    ref2 = a2d[i];
    while (i < j)
    {
	while (i < j && ref -  alf[j] < 0)
	    j--;
	if (i != j) {
	   a2d[i] = a2d[j];
	    alf[i++] = alf[j];
	}
	while (i < j && ref - alf[i] > 0)
	    i++;
	if (i != j) {
	   a2d[j] = a2d[i];
	    alf[j--] = alf[i];
	}
    }
    a2d[j] = ref2;
    alf[j] = ref;
    if (left < --j)
	quicksort2lfd(left, j);
    if (++i < right)
	quicksort2lfd(i, right);
}
/**********************************************************************/

/* sort array into INCREASING ORDER, and shuffle array2 at the same time
   if array2 = 0,1,...,N and A represents the original _unsorted_ array,
   then after the sort, A[array2[0]], A[array2[1]], ...
   represents the sorted data in array.
*/
void sort2lfd(int num_elements, double *array, int *array2)
{
    if (num_elements < 2)
	return;
    alf = array;
    a2d = array2;
    quicksort2lfd(0, num_elements - 1);
}

/**********************************************************************/
static void quicksort1lf(int left, int right)
{
    long i, j;
    double ref;
    int ref2;
    i = left;
    j = right;
    ref = alf[i];
    while (i < j)
    {
	while (i < j && ref -  alf[j] < 0)
	    j--;
	if (i != j) {
	    alf[i++] = alf[j];
	}
	while (i < j && ref - alf[i] > 0)
	    i++;
	if (i != j) {
	    alf[j--] = alf[i];
	}
    }
    alf[j] = ref;
    if (left < --j)
	quicksort1lf(left, j);
    if (++i < right)
	quicksort1lf(i, right);
}
/**********************************************************************/

void sort1lf(int num_elements, double *array)
{
    if (num_elements < 2)
	return;
    alf = array;
    quicksort1lf(0, num_elements - 1);
}

/**********************************************************************/

static void quicksort2fd(int left, int right)
{
    long i, j;
    float ref;
    int ref2;
    i = left;
    j = right;
    ref = af[i];
    ref2 = a2d[i];
    while (i < j)
    {
	while (i < j && ref -  af[j] < 0)
	    j--;
	if (i != j) {
	   a2d[i] = a2d[j];
	    af[i++] = af[j];
	}
	while (i < j && ref - af[i] > 0)
	    i++;
	if (i != j) {
	   a2d[j] = a2d[i];
	    af[j--] = af[i];
	}
    }
    a2d[j] = ref2;
    af[j] = ref;
    if (left < --j)
	quicksort2fd(left, j);
    if (++i < right)
	quicksort2fd(i, right);
}
/**********************************************************************/

/* sort array into INCREASING ORDER, and shuffle array2 at the same time
   if array2 = 0,1,...,N and A represents the original _unsorted_ array,
   then after the sort, A[array2[0]], A[array2[1]], ...
   represents the sorted data in array.
*/
void sort2fd(int num_elements, float *array, int *array2)
{
    if (num_elements < 2)
	return;
    af = array;
    a2d = array2;
    quicksort2fd(0, num_elements - 1);
}

/**********************************************************************/
static void quicksort1f(int left, int right)
{
    long i, j;
    float ref;
    int ref2;
    i = left;
    j = right;
    ref = af[i];
    while (i < j)
    {
	while (i < j && ref -  af[j] < 0)
	    j--;
	if (i != j) {
	    af[i++] = af[j];
	}
	while (i < j && ref - af[i] > 0)
	    i++;
	if (i != j) {
	    af[j--] = af[i];
	}
    }
    af[j] = ref;
    if (left < --j)
	quicksort1f(left, j);
    if (++i < right)
	quicksort1f(i, right);
}
/**********************************************************************/

void sort1f(int num_elements, float *array)
{
    if (num_elements < 2)
	return;
    af = array;
    quicksort1f(0, num_elements - 1);
}

/*
Local Variables:
compile-command: "g++ -c -g utiltkm.cc"
End:
*/

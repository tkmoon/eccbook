// An interface to quicksort

/**********************************************************************/
static double *adouble;
static int *aint;
static int *a2d;
static float *af;


/**********************************************************************/

static void quicksort1double(int left, int right)
{
   long i, j;
   double ref;
   int ref2;
   i = left;
   j = right;
   ref = adouble[i];
   while (i < j) {
	  while (i < j && ref -  adouble[j] < 0)
		 j--;
	  if (i != j) {
		 adouble[i++] = adouble[j];
	  }
	  while (i < j && ref - adouble[i] > 0)
		 i++;
	  if (i != j) {
		 adouble[j--] = adouble[i];
	  }
   }
   adouble[j] = ref;
   if (left < --j)
	  quicksort1double(left, j);
   if (++i < right)
	  quicksort1double(i, right);
}

void sort1double(int num_elements, double *array)
{
    if (num_elements < 2)
	return;
    adouble = array;
    quicksort1double(0, num_elements - 1);
}


static void quicksort2double_int(int left, int right)
{
   long i, j;
   double ref;
   int ref2;
   i = left;
   j = right;
   ref = adouble[i];
   ref2 = a2d[i];
   while (i < j) {
	  while (i < j && ref -  adouble[j] < 0)
		 j--;
	  if (i != j) {
		 a2d[i] = a2d[j];
		 adouble[i++] = adouble[j];
	  }
	  while (i < j && ref - adouble[i] > 0)
		 i++;
	  if (i != j) {
		 a2d[j] = a2d[i];
		 adouble[j--] = adouble[i];
	  }
   }
   a2d[j] = ref2;
   adouble[j] = ref;
   if (left < --j)
	  quicksort2double_int(left, j);
   if (++i < right)
	  quicksort2double_int(i, right);
}

static void quicksort2int_int(int left, int right)
{
   long i, j;
   double ref;
   int ref2;
   i = left;
   j = right;
   ref = aint[i];
   ref2 = a2d[i];
   while (i < j) {
	  while (i < j && ref -  aint[j] < 0)
		 j--;
	  if (i != j) {
		 a2d[i] = a2d[j];
		 aint[i++] = aint[j];
	  }
	  while (i < j && ref - aint[i] > 0)
		 i++;
	  if (i != j) {
		 a2d[j] = a2d[i];
		 aint[j--] = aint[i];
	  }
   }
   a2d[j] = ref2;
   aint[j] = ref;
   if (left < --j)
	  quicksort2int_int(left, j);
   if (++i < right)
	  quicksort2int_int(i, right);
}

void sort2double_int(int num_elements, double *array, int *array2)
{
    if (num_elements < 2)
	return;
    adouble = array;
    a2d = array2;
    quicksort2double_int(0, num_elements - 1);
}

void sort2int_int(int num_elements, int *array, int *array2)
{
    if (num_elements < 2)
	return;
    aint = array;
    a2d = array2;
    quicksort2int_int(0, num_elements - 1);
}

/*
Local Variables:
compile-command: "clang++ -c quicksortstuff.cc"
End:
*/

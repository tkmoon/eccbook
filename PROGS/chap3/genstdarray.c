/**************************************
*
*  Program: genstdarray -- generate a standard array for a code
*  (used for example purposes)
*  This is not very general, having aspects hard-coded, but
*  suffices to generate some examples.
*
*  It also produces the weight distribution.
*
*  This can be readily modified for other codes, but 
*  since the computations are all explicit, this is very 
*  impractical for large n or k.
*
*  Todd K. Moon
*  Utah State University
*
*  Date: 1/8/03
*
***************************************/


#define N 7
#define K 3

void n2bin(int m, char *mv, int k)  /* convert a number to binary */
{
  int i, mask;
  for(i = 0,mask=1; i < k; i++, mask<<=1) {
	mv[i] = (m & mask) != 0;
  }
}

void printbin(int e1, int n)  /* print a binary number */
{
  int i, mask;
  for(i = 0,mask=1; i < n; i++, mask<<=1) {
	if(e1 & mask) printf("1"); else printf("0");
  }
}

void updatewt(int *wtdist, char *cv, int n)
{
  int wt=0, i;
  for(i = 0; i < n; i++) {
	wt += cv[i];
  }
  wtdist[wt]++;
}

void encode(char *mv, char *cv, char G[][K], int n, int k) /* encode using G */
{
  int i,j;
  char sum;
  for(i = 0; i < n; i++) {
	sum = 0;
	for(j = 0; j < k; j++) {
	  sum += G[i][j]*mv[j];
	}
	cv[i] = sum % 2;
  }
}

void encode2(char *ev, char *sv, char H[][N], int nk, int n) /* compute syndrome */
{
  int i,j;
  char sum;
  for(i = 0; i < nk; i++) {
	sum = 0;
	for(j = 0; j < n; j++) {
	  sum += H[i][j]*ev[j];
	}
	sv[i] = sum % 2;
  }
}

void bsc(char *cv,int e1,int *r,int n)     /* add on row leader (error) */
{
  int i,mask;
  *r = 0;
//  for(i = 0,mask=(1<<(n-1)); i < n; i++,mask >>= 1) {
  for(i = 0,mask=1; i < n; i++,mask <<= 1) {
	*r |= cv[i]*mask;
	if(e1 & mask) *r ^= mask;
  }
}

void printvec(char *cv, int n)
{
  int i;
  for(i = 0; i < n; i++) printf("%d",cv[i]);
}

int findminwt(char *used,int n) /* find the unused n-tuple of lightest weight */
{
  int i,minwt,mini;
  int j,wt;

  minwt = n;
  for(i = 0; i < (1<<n); i++) {
	if(used[i]) continue;
	for(wt=0, j = 0; j < n; j++) { /* compute the weight of the number i */
	  if(i & (1<<j)) wt++;
	}
	if(wt < minwt) {
	  minwt = wt;
	  mini = i;
	}
  }
  return mini;
}
	


main() 
{
  int n = N;
  int k = K;
  char mv[K];  /* message vector */
  char cv[N];  /* code vector */
  char sv[N-K];/* syndrome vector */

  char G[N][K] = {{0,1,1},
				  {1,0,1},
			      {1,1,0},
				  {1,1,1},
				  {1,0,0},
				  {0,1,0},
                  {0,0,1}};

  char H[N-K][N] = {{1,0,0,0,  0,1,1}, 
					{0,1,0,0,  1,0,1}, 
                    {0,0,1,0,  1,1,0}, 
                    {0,0,0,1,  1,1,1}};
  int e[(1<<(N-K))];			/* save the coset leaders */
  char used[1<<N];				/* indicate which are used */
  int wtdist[N];				/* weight distribution array */

  int twok, twon;
  int i,j;
  int e1,m,y;
  char ev[N];
  int first = 1;

  twon = 1<<n;   /* 2^n */
  twok = 1<<k;   /* 2^k */

  for(i = 0; i < twon; i++) used[i] = 0;
  for(i = 0; i < N; i++) wtdist[i] = 0;

  for(i = 0; i < (1<<N-K); i++) { /* for each row */
	if(!first) {				/* if not first time */
	  e1 = findminwt(used,N);	/* find column leader of lightest weight */
	}
	else {						/* first time through, column leader=0 */
	  e1 = 0;
	}
	e[i] = e1;
	for(m = 0; m < twok; m++) {  /* over all input vectors */
	  n2bin(m,mv,k);
	  encode(mv,cv,G,n,k);		/* code the vector */
	  if(first) updatewt(wtdist,cv,n);	/* compute the weight distribution */
	  bsc(cv,e1,&y,n);			/* add on row leader (error) */
	  used[y] = 1;				/* set that this one is used */
	  printbin(y,n);
	  if(m < twok-1) printf("&"); else {printf("\\\\"); printf("\n"); }
	}
	first = 0;
  }
  printf("\n\n\n");
  for(i = 0; i < (1<<(N-K)); i++) {
	e1 = e[i];
    printbin(e1,n); printf("& ");
	n2bin(e1,ev,n);
	encode2(ev,sv,H,n-k,n);
	printvec(sv,n-k); printf("\\\\ \n");
  }

  /* print out the weight distribution */
  printf("\n\nWeight distribution:\n");
  for(i = 0; i < N; i++) {
	printf("(%d,%d)\n",i,wtdist[i]);
  }

}

/*
Local Variables:
compile-command:"gcc -o genstdarray -g genstdarray.c"
End:
*/

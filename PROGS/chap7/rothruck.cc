// rothruck.cc -- Given a polynomial in two variables Q(x,y),
// determine all the y roots up to degree D using the Roth-Ruck 
// algorithm

// Todd K. Moon, February 11, 2004

#ifdef GFTYPE
#include "GFNUM2m.h"
#endif
#include "rothruck.h"

static rpolynode *rpolystart = NULL;   // pointer to start of list
static rpolynode *rpolylist = NULL; // pointer to list as built
static int rothD;				// maximum degree sought
static int rothnodenum; 		// counts the nodes in tree
static int numinrrpool = 0;		// number of solutions found

#define DO_PRINT


void
QxytoQxxya(const polynomialT<polynomialT<TYPE> > &Q, const TYPE &a,
		   polynomialT<polynomialT<TYPE> > &P);
void rothrucktree(polynomialT<polynomialT<TYPE> >  &Qu, int u, 
					 TYPE *f);
void RootList(TYPE *R,const polynomialT<TYPE> &Q,int &numroots);
int binom(int n, int k);



rpolynode *
rothruck(polynomialT<polynomialT<TYPE> >&Q,int D)
{


   int u = 0;
   TYPE *f = new TYPE[D+1];		// space for the output polynomial

   // build the list where the answers will go
   numinrrpool = 0;
   rpolylist = rpolystart;		// point to beginning of list

   rothD = D;					// set the maximum degree
   rothnodenum = 0;

   // call the function which does the recursive work
   rothrucktree(Q,u,f);

   // Dump the output
   // cout << "All the factors:\n";
//    rpolynode *tptr;
//    for(tptr = rpolystart; tptr != NULL; tptr = tptr->next) {
// 	  cout << tptr->f << endl;
//    }
   return rpolystart;
}

void 
rothrucktree(polynomialT<polynomialT<TYPE> > &Qu, int u, TYPE *f)
{

   rothnodenum++;
#ifdef DO_PRINT
   cout << "nodenum=" << rothnodenum << endl;

   cout <<"Qu=" << Qu << endl;
#endif
   int i;
   if(Qu[0] == 0) {
#ifdef DO_PRINT
	  cout << "***** f=";
	  for(i = 0; i < u; i++) {
		 cout << f[i] << " ";
	  }
	  cout << " ******" << endl;
#endif
	  // add the polynomial to the list
	  if(rpolystart == NULL) {	// get the list started
		 rpolystart = new rpolynode[1];
		 rpolylist = rpolystart;
		 rpolylist->next = NULL;
	  }
	  else if(rpolylist->next == NULL) { // need to add to list
		 rpolylist->next = new rpolynode[1];
		 rpolylist = rpolylist->next;
		 rpolylist->next = NULL;
	  }
	  (rpolylist->f).setc(f,u-1);
   }
   else if(u <= rothD) {  // try another branch
	  int qdeg = Qu.getdegree();
	  TYPE *R = new TYPE[qdeg];  // possibly this many roots
	  int numroots;
	  // copy Q(0,y) into a single-variable polynomial for speed
	  polynomialT<TYPE> *Q0y = polytemp<TYPE>::gettemppoly(qdeg);
	  for(i = 0; i <= qdeg; i++) {
		 Q0y->operator[](i) = Qu[i][0];
	  }
	  RootList(R,*Q0y,numroots);

#ifdef DO_PRINT
	  cout << "tree depth=" << u << endl;
	  cout << "Q(0,y)=" << *Q0y << endl;
	  cout <<"Roots: ";
	  for(i = 0; i < numroots; i++){
		 cout << R[i] << " ";
	  }
	  cout << endl;
#endif
	  polynomialT<polynomialT<TYPE> > Qv;
	  for(i = 0; i < numroots; i++) {
#ifdef DO_PRINT
		 cout << "Root=" << R[i] << endl;
#endif
		 QxytoQxxya(Qu,  R[i], Qv);
		 f[u] = R[i];
		 rothrucktree(Qv, u+1, f);
	  }
	  delete [] R;
   }
   else {
	  cout << "No output\n";
   }
}


void RootList(TYPE *R,const polynomialT<TYPE> &Q,int &numroots)
{
   // this is a rather primitive approach to finding all the 
   // roots.  A better approach would be a Chien search, but
   // this serves the purposes of demonstration.
   int i = 0;
   numroots = 0;
#ifdef MODARTYPE
   int m = Q[0].getm();			// get the size of the field
   for(i = 0; i < m; i++) {
	  if(Q(i) == 0) {
		 R[numroots++] = i;
	  }
   }
#endif
#ifdef GFTYPE
   if(Q(0) == 0) R[numroots++] = 0;
   for(i = 0; i < Q[0].getN(); i++) {
	  if(Q(A^i) == 0) R[numroots++] = A^i;
   }
#endif  
}


void
QxytoQxxya(const polynomialT<polynomialT<TYPE> > &Q, const TYPE &a,
		   polynomialT<polynomialT<TYPE> > &P)
{
   int i,k,j;
   TYPE nk;

   // This is a rather literal translation.  Probably more efficiency
   // could be obtained if the polynomials were represented in
   // a kind of matrix form.  But this serves for purposes of
   // demonstration.
   int mod = a.character();
   int ydeg = Q.getdegree();
   P = Q;						// to get started
   for(i = 1; i <= ydeg; i++) P[i] = polynomialT<TYPE>(0);
   //cout << "i=0  P=" << P << endl;
   for(i = 1; i <= ydeg; i++) {
	  for(k = 0; k <= i; k++) {
		 nk = TYPE(binom(i,k) % mod);
		 if(nk != 0) {
			P[k] += ((Q[i] << k)*(a^(i-k)))*nk;
		 }
	  }
	  // cout <<"i=" << i << "   P=" << P << endl;
   }
   // now find the largest m such that x^m | P
   int m = P[0].getdegree();		// starting point
   for(i = 0; i <= ydeg; i++) {
	  for(j = 0; j <= P[i].getdegree(); j++) {
		 if(P[i][j] != 0){		// nonzero coefficient found
			if(j < m) m = j;
			break;
		 }
	  }
   }
   // cout << "P=" << P << endl;
   // cout << "m=" << m << endl;
   // now divide by x^m
   for(i = 0; i <= ydeg; i++) {
	  P[i] >>= m;
   }
   // cout << "P/x^m=" << P << endl;
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

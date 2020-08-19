// Compute kotters interpolation algorithm.
// Given points (xi,yi) for i=0,...,n-1, a weighted degree in
// wlex, and an interpolation multiplicity mi, return the polynomial
// Q(x,y) interpolating the points.

// If mi==NULL, then a single order m1 is used at each point.

// This implementation is a rather slavish adherance to the 
// statement of the algorithm, using bi-variate polynomials.
// It could no doubt be made much more efficient by
// - computing the w-degree more efficiently
// - computing the Hasse derivatives more efficiently
// - implementing less literally the polynomial arithmetic

// Todd K. Moon, February 12, 2004

#define DO_PRINT

#define BIGINT 65535


polynomialT<polynomialT<TYPE> >
kotter(int n,int L,TYPE *xi, TYPE *yi, int *mi,int *wdeg,int rlex,
	   int m1)
   // if mi=NULL, then the same multiplicity is used by every point
{


   polynomialT<polynomialT<TYPE> > f;
   polynomialT<polynomialT<TYPE> > gj[L+1];
   TYPE lc[] = {1,1};
   polynomialT<TYPE> l1(1,lc); // x + 1 --- place holder to compute

 
   int computewdeg(const polynomialT<polynomialT<TYPE> > &Q,const int *wdeg,
				   int &xi, int &yj);
   TYPE computeD(int r,int s,const polynomialT<polynomialT<TYPE> > &Q,
					TYPE a, TYPE b);
   void findminwdeg(const polynomialT<polynomialT<TYPE> > *gj,
				 int *Jlist,int numinJ,int &jstar,int &xistar,int &yjstar,
				 int *wdeg,int rlex);


   int xi1,yj1;
   int i,j,j1,jstar,xistar,yjstar;
   polynomialT<TYPE> Deltaj[L+1];  // set these as polynomials,
   polynomialT<TYPE> Delta;        // to avoid having to cast later
   int Jlist[L+1];
   int numinJ,wd,wdmin;
   int r, s;
   int rs;
   int m;

   // Initialize
   for(int j = 0; j <= L; j++) {
	  gj[j].setc(j,polynomialT<TYPE>(1)); // make a y^j
	  gj[j].setvarname("y");
	  gj[j].setbeforeafterprint("y","(",")");
   }
   cout << endl;
   // compute the order of the interpolator
   int C = 0;
   for(i = 0; i < n; i++) {	 // loop through the points to interpolate

	  if(mi==NULL) {
		 m=m1;
	  }
	  else {
		 m = mi[i];
	  }
	  C = (m+1)*m/2; // number of interpolation conditions
#ifdef DO_PRINT
	  cout << "i=" << i << "  C=" << C << endl;
#endif
	  for(rs= 0; rs < C; rs++) { // walk through in lex order
		 s = rs % m;
		 r = rs / m;
#ifdef DO_PRINT
		 cout << "   (r,s)=" << r << "," << s << endl;
	  for(j = 0; j <= L; j++) {
		 cout << "   gj[" << j << "]=" << gj[j] << endl;
	  }
#endif
		 numinJ = 0;			// start from a fresh J
		 for(j = 0; j <= L; j++) {
			// compute jth discrepancy
			Deltaj[j][0] = computeD(r,s,gj[j],xi[i],yi[i]);
			if(Deltaj[j][0] != 0) {
			   Jlist[numinJ++] = j;
			}
		 }
		 if(numinJ) {			// J != emptyset
			findminwdeg(gj,Jlist,numinJ,jstar,xistar,yjstar,wdeg,rlex);
			f = gj[jstar];
			Delta[0] = Deltaj[jstar][0];
#ifdef DO_PRINT
			cout << "      J= ";
			for(j1 = 0; j1 < numinJ; j1++)cout<< Jlist[j1]<<" "; cout << endl;
			cout << "      ";
			for(j1 = 0; j1 < numinJ; j1++)
			   cout << "Delta[" << Jlist[j1] << "]=" <<
				  Deltaj[Jlist[j1]][0] << " ";
			cout << endl;
			cout << "      jstar= " << jstar << "  f=" << f << "  Delta=" <<
			   Delta[0] << endl;
#endif
			for(j1 = 0; j1 < numinJ; j1++) { // for j in J
			   j = Jlist[j1];
			   if(j != jstar) {
				  gj[j] = gj[j]*Delta - f*Deltaj[j];
			   }
			   else if(j == jstar) {
				  l1[0] = -xi[i];
				  gj[j] = f*l1;
			   }
			   gj[j].setvarname("y");
			}  // for j in J
		 }  // if J != emptyset
	  } // for rs
   } // for i

   // find the one of minimum weight
#ifdef DO_PRINT
   cout << "Final polynomials:" << endl;
	  for(j = 0; j <= L; j++) {
		 int xi, yj;
		 cout << "   gj[" << j << "]=" << gj[j] << "  wdeg=" <<
			computewdeg(gj[j], wdeg, xi,yj) << endl;
	  }
#endif
   for(j = 0; j <= L; j++) Jlist[j] = j;
   numinJ = L+1;
   findminwdeg(gj,Jlist,numinJ,jstar,xistar,yjstar,wdeg,rlex);
   return gj[jstar];
}

TYPE evaluate(polynomialT<polynomialT<TYPE> > &Q,TYPE a,
				 TYPE b)
// Find Q(a,b)
{
   TYPE e(0);				// evaluation value
   polynomialT<TYPE> x;
   for(int i = 0; i <= Q.getdegree(); i++) {
	  x = Q[i];
	  e += x(a) * (b^i);
   }
   return e;
}



void findminwdeg(const polynomialT<polynomialT<TYPE> > *gj,
				 int *Jlist,int numinJ,int &jstar,int &xistar,int &yjstar,
				 int *wdeg,int rlex)
{
   int wdmin,j1,j,wd,xi,yj;
   int computewdeg(const polynomialT<polynomialT<TYPE> > &Q,const int *wdeg,
				   int &xi, int &yj);

   wdmin = BIGINT;
   for(j1 = 0; j1 < numinJ; j1++) {
	  j = Jlist[j1];
	  wd = computewdeg(gj[j],wdeg,xi,yj);
	  if(wd < wdmin) {
		 jstar = j;
		 xistar = xi;
		 yjstar = yj;
		 wdmin = wd;
	  }
	  else if(wd == wdmin) { // if the same, go on the basis of lex
		 if(rlex==0) {	// lex order
			if(xi<xistar) { // new one is lower
			   jstar = Jlist[j];
			   xistar = xi;
			   yjstar = yj;
			}
		 }
		 else {		// rlex order
			if(xi > xistar) { // new one is lower
			   jstar = Jlist[j];
			   xistar = xi;
			   yjstar = yj;
			}
		 }
	  }
   }
}



int computewdeg(const polynomialT<polynomialT<TYPE> > &Q,const int *wdeg,
				int &xi, int &yj)
{
   // this is an inefficient search, being exhaustive.
   // However, it gets the job done.
   int degy = Q.getdegree();
   int degx;
   int deg;
   int maxdeg = 0;
   xi = yj = 0;
   for(int i = 0; i <= degy; i++) {
	  degx = Q[i].getdegree();
	  for(int j = 0; j <= degx; j++) {
		 if(Q[i][j] != 0) {
			deg = wdeg[0]*j + wdeg[1]*i;
			if(deg > maxdeg) {
			   maxdeg = deg;
			   xi = j;
			   yj = i;
			}
		 }
	  }
   }
   return maxdeg;
}


TYPE computeD(int r,int s,const polynomialT<polynomialT<TYPE> > &Q,
				 TYPE a, TYPE b)
{
   int degy = Q.getdegree();
   TYPE Drs;
   int degx;
   int jchs, ichr;
   int binommod(int n, int k, int m);
   TYPE aij;
   int charac = a.character();

   for(int j = s; j <= degy; j++) {
	  jchs = binommod(j,s,charac);
	  if(jchs) {
		 degx = Q[j].getdegree();
		 for(int i = r; i <= degx; i++) {
			aij = Q[j][i];
			if(aij != 0) {
			   ichr = binommod(i,r,charac);
			   if(ichr) {
				  Drs += aij*(a^(i-r))*(b^(j-s));
			   }
			}
		 }
	  }
   }
   return Drs;
}

int binommod(int n, int k, int character)
{
   double prod = 1;
   long int ipod;
   double i,j;
   if(k > n) return 0;
   if(k > n/2) k = n-k;
   if(k <= 1) {
	  if(k == 0) return 1;
	  if(k == 1) return n % character;
   }
   for(i = n-k+1, j=1; i <= n; i++,j++) {
	  prod *= i/j;
   }
   ipod = (long int)(prod+0.5);
   return ipod % character;
}

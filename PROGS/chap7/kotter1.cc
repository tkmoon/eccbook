// Compute kotters interpolation algorithm.
// Given points (xi,yi) for i=0,...,n-1, a weighted degree in
// wlex.  This verion uses m=1, then returns the codeword

// Todd K. Moon, February 12, 2004

#define DO_PRINT

#define BIGINT 65535

TYPE evaluate(polynomialT<polynomialT<TYPE> > &Q,TYPE a,TYPE b);

polynomialT<TYPE>
kotter1(int n,int k,TYPE *xi, TYPE *yi, int *wdeg,int rlex)
{
   polynomialT<polynomialT<TYPE> > f;
   polynomialT<polynomialT<TYPE> > gj[2];
   TYPE lc[] = {1,1};
   polynomialT<TYPE> l1(1,lc); // x + 1 --- place holder to compute

 
   int computewdeg(const polynomialT<polynomialT<TYPE> > &Q,const int *wdeg,
				   int &xi, int &yj);
   TYPE computeD(int r,int s,const polynomialT<polynomialT<TYPE> > &Q,
					TYPE a, TYPE b);
   void findminwdeg(const polynomialT<polynomialT<TYPE> > *gj,
				 int *Jlist,int numinJ,int &jstar,int &xistar,int &yjstar,
				 int *wdeg,int rlex);

   int L = 1;
   int xi1,yj1;
   int i,j,j1,jstar,xistar,yjstar;
   polynomialT<TYPE> Deltaj[L+1];  // set these as polynomials,
   polynomialT<TYPE> Delta;        // to avoid having to cast later
   int Jlist[L+1];
   int numinJ,wd,wdmin;

   // Initialize
   for(int j = 0; j <= L; j++) {
	  gj[j].setc(j,polynomialT<TYPE>(1)); // make a y^j
	  gj[j].setvarname("y");
	  gj[j].setbeforeafterprint("y","(",")");
   }
   cout << endl;
   for(i = 0; i < n; i++) {	 // loop through the points to interpolate
#ifdef DO_PRINT
	  cout << "i=" << i << endl;
#endif
#ifdef DO_PRINT
	  for(j = 0; j <= L; j++) {
		 cout << "   gj[" << j << "]=" << gj[j] << endl;
	  }
#endif
	  numinJ = 0;			// start from a fresh J
	  for(j = 0; j <= L; j++) {
		 // compute jth discrepancy
		 Deltaj[j][0] = evaluate(gj[j],xi[i],yi[i]);
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
   polynomialT<TYPE> P1(gj[jstar][1]);
   polynomialT<TYPE> P0(-gj[jstar][0]);
#ifdef DO_PRINT
   cout << "P0=" << P0 << endl;
   cout << "P1=" << P1 << endl;
#endif
   polynomialT<TYPE> r = P0 % P1;
   polynomialT<TYPE> p;
   if(r == 0) {  
	  p = polynomialT<TYPE>::getlastquotient();
	  if(p.getdegree() > k-1) {
		 cout << "Uncorrectable error pattern\n";
		 p = TYPE(0);
	  }
   }
   else {
	  cout << "Uncorrectable error pattern\n";
	  p = TYPE(0);
   }
   return p;
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

/*
Local Variables:
compile-command: "g++ -c -g kotter1.cc"
End:
*/

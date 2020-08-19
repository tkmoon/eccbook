
#include <iostream>
using namespace std;

static int doprint = 1;
static int printmatlabform = 0;

static int getdoprintvalue(void)
{
   return doprint;
}

static void setdoprintvalue(int printvalue)
{
   doprint = printvalue;
}

static void setprintmatlabform(int matlabform)
{
   printmatlabform = matlabform;
}


static void printdoublevec(string s, double *vec, int N)
{
   if(doprint) {
	  if(!printmatlabform) cout << s << " ";
	  else cout << s << " = [";
	  for(int i = 0; i < N; i++)  cout << vec[i] << " ";
	  if(printmatlabform) cout << "];";
	  cout << endl;
   }
}


static void printintmat(string s, int **mat, int m, int n)
{
   if(!doprint) return;
   if(!printmatlabform) cout << s << " ";
   else cout << s << " = [";
   cout << endl;
   for(int i = 0; i < m; i++) {
	  for(int j = 0; j < n; j++) {
		 cout << mat[i][j] << " ";
	  }
	  if(printmatlabform) {
		 if(i == m-1) {
			cout << "];" << endl;
		 }
		 else {
			cout << endl;
		 }
	  }
	  else {
		 cout << endl;
	  }
   }
}


static void printdoublemat(string s, double **mat, int m, int n)
{
   if(!doprint) return;
   if(!printmatlabform) cout << s << " ";
   else cout << s << " = [";
   cout << endl;
   for(int i = 0; i < m; i++) {
	  for(int j = 0; j < n; j++) {
		 cout << mat[i][j] << " ";
	  }
	  if(printmatlabform) {
		 if(i == m-1) {
			cout << "];" << endl;
		 }
		 else {
			cout << endl;
		 }
	  }
	  else {
		 cout << endl;
	  }
   }
}

static void printsignedcharmat(string s, signed char **mat, int m, int n)
{
   if(!doprint) return;
   if(!printmatlabform) cout << s << " ";
   else cout << s << " = [";
   cout << endl;
   for(int i = 0; i < m; i++) {
	  for(int j = 0; j < n; j++) {
		 cout << int(mat[i][j]) << " ";
	  }
	  if(printmatlabform) {
		 if(i == m-1) {
			cout << "];" << endl;
		 }
		 else {
			cout << endl;
		 }
	  }
	  else {
		 cout << endl;
	  }
   }
}

static void printunsignedcharmat(string s, unsigned char **mat, int m, int n)
{
   if(!doprint) return;
   if(!printmatlabform) cout << s << " ";
   else cout << s << " = [";
   cout << endl;
   for(int i = 0; i < m; i++) {
	  for(int j = 0; j < n; j++) {
		 cout << int(mat[i][j]) << " ";
	  }
	  if(printmatlabform) {
		 if(i == m-1) {
			cout << "];" << endl;
		 }
		 else {
			cout << endl;
		 }
	  }
	  else {
		 cout << endl;
	  }
   }
}

static void printboolvec(string s, bool *vec, int m)
{
   if(doprint) {
	  if(!printmatlabform) cout << s << " ";
	  else cout << s << " = [";
	  for(int i = 0; i < m; i++)  cout << int(vec[i]) << " ";
	  if(printmatlabform) cout << "];";
	  cout << endl;
   }

}

static void printboolmat(string s, bool **mat, int m, int n)
{
   if(!doprint) return;
   if(!printmatlabform) cout << s << " ";
   else cout << s << " = [";
   cout << endl;
   for(int i = 0; i < m; i++) {
	  for(int j = 0; j < n; j++) {
		 cout << int(mat[i][j]) << " ";
	  }
	  if(printmatlabform) {
		 if(i == m-1) {
			cout << "];" << endl;
		 }
		 else {
			cout << endl;
		 }
	  }
	  else {
		 cout << endl;
	  }
   }
}

static void printbitmat(string s, BITTYPE **mat, int m, int n)
{
   if(!doprint) return;
   if(!printmatlabform) cout << s << " ";
   else cout << s << " = [";
   cout << endl;
   for(int i = 0; i < m; i++) {
	  for(int j = 0; j < n; j++) {
		 cout << int(mat[i][j]) << " ";
	  }
	  if(printmatlabform) {
		 if(i == m-1) {
			cout << "];" << endl;
		 }
		 else {
			cout << endl;
		 }
	  }
	  else {
		 cout << endl;
	  }
   }
}

static void printintvec(string s, int *vec, int N)
{
   if(doprint) {
	  if(!printmatlabform) cout << s << " ";
	  else cout << s << " = [";
	  for(int i = 0; i < N; i++)  cout << vec[i] << " ";
	  if(printmatlabform) cout << "];";
	  cout << endl;
   }
}


static void printsignedcharvec(string s, signed char *vec, int N)
{
   if(doprint) {
	  if(!printmatlabform) cout << s << " ";
	  else cout << s << " = [";
	  for(int i = 0; i < N; i++)  cout << int(vec[i]) << " ";
	  if(printmatlabform) cout << "];";
	  cout << endl;
   }
}

static void printbitvec(string s, BITTYPE *vec, int N)
{
   if(doprint) {
	  if(!printmatlabform) cout << s << " ";
	  else cout << s << " = [";
	  for(int i = 0; i < N; i++)  cout << int(vec[i]) << " ";
	  if(printmatlabform) cout << "];";
	  cout << endl;
   }
}




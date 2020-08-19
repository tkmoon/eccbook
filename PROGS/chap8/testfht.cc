//
//
//  Program: testfht.cc --- test the fast Hadamard transform routine
//
//
//  Todd K. Moon
//  Utah State University
//
//  Date:  March 24, 2004

#include <iostream>
using namespace std;
void fht(int *F, int m);
void fhtinv(int *F, int m);

int main()
{
   int m = 3;
   int n = (1<<m);
   int F[] = {-1,1,1,-1,1,1,-1,1};
   int i;
   cout << "Original:\n";
   for(i = 0; i < n; i++) {
	  cout << F[i] << " ";
   }
   cout << endl;
   

   fht(F,3);
   cout << "transform:\n";
   for(i = 0; i < n; i++) {
	  cout << F[i] << " ";
   }
   fhtinv(F,3);
   cout << "inverse transform:\n";
   for(i = 0; i < n; i++) {
	  cout << F[i] << " ";
   }
   
   cout << endl;
   cout << endl;
}



/*
Local Variables:
compile-command: "g++ -o testfht testfht.cc fht.cc"
End:
*/



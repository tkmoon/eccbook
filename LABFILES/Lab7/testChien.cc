//  Program: testChienSearch.cc
//  Todd K. Moon

#include "ChienSearch.h"



main()
{
   GFNUM2m::initgf(4,0x13);  // 1 011 = 1+x+x^4
   ChienSearch Searcher(3);
   GFNUM2m Lambda[] = {1,1,0,A^5}; // 1+x + A^5 x^3
   GFNUM2m *p = Searcher.Search(Lambda,3);
   for(int i = 0; i < Searcher.getnroots(); i++) {
	  cout << p[i] << " ";
   }
   cout << endl;
}
   

/*
Local Variables:
compile-command: "g++ -o testChien -g testChien.cc ChienSearch.cc GFNUM2m.cc"
End:
*/



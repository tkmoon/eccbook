//  Program: maketable.cc
//  Make the addition/multiplication table

#include "GFNUM2m.h"

int main()
{
   int g = 0x13;  // 1 0011 = 1+d+d^4
   GFNUM2m::initgf(4,g);

//    int g = 0xB;  // 1011 = 1+d+d^3
//    GFNUM2m::initgf(3,g);

//	int g = 0xD;   // 1101 = 1+d^2+d^3
//    GFNUM2m::initgf(3,g);

   GFNUM2m a,b,c;

   int i,j;

   cout << "log: ";
   for(i = 0; i <= GFNUM2m::getN(); i++) {
	  a = i;
	  cout << a.getp() << " ";
   }
   cout << endl;
   cout << "zlog: ";
   for(i = 1; i <= GFNUM2m::getN(); i++) {
	  a = i;
	  b = 1+a;
	  // b = 1+(A^i);
	  cout << " " << b.getp();
   }
   cout << endl;
   for(i = 0; i <= GFNUM2m::getN(); i++) {
	  a = i;
	  cout << a << " " << a.getv() << "    ";
	  for(j = 0; j <= GFNUM2m::getN(); j++) {
		 b = j;
		 if(j < i) {			// addition table
			c = a+b;
			cout << c.getv();
		 }
		 else {					// multiplication table
			c = a*b;
			cout << c.getv();
		 }
		 cout << " ";
		 if(j==i-1) cout << "|";
	  }
	  cout << endl;
   }


   // now the addition table in terms of alpha
   cout << "\n\n\n\n";
   cout << "&";
   for(i = 0; i <= GFNUM2m::getN(); i++) {
	  if(i==0) a = 0;
	  else a = A^(i-1);
	  if(a==0) cout << "0";
	  else if(a==1) cout << "1";
	  else if(a==A) cout << "\\alpha";
	  else
		 cout << "\\alpha^{" << a.getp() << "}";
	  cout  << "& ";
   }
   cout << endl;
   for(i = 0; i <= GFNUM2m::getN(); i++) {
	  if(i==0) a = 0;
	  else a = A^(i-1);

	  if(a==0) cout << "0";
	  else if(a==1) cout << "1";
	  else if(a==A) cout << "\\alpha";
	  else
		 cout << "\\alpha^{" << a.getp() << "}";
	  cout  << "& ";

	  for(j = 0; j <= GFNUM2m::getN(); j++) {
		 if(j==0) b = 0;
		 else b = A^(j-1);
		 c = a+b;
		 if(c ==0) cout << "0";
 		 else if(c==1) cout << "1";
		 else if(c==A) cout << "\\alpha";
 		 else
 			cout << "\\alpha^{" << c.getp() << "}";
		 cout  << "& ";
	  }
	  cout << "\n";
   }

}

/*
Local Variables:
compile-command: "g++ -g -o maketable maketable.cc GFNUM2m.cc"
End:
*/



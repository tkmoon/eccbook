/**************************************
*
*  polarstack.h -- a simple stack non-templatized
*      stack for the polar decoder
*
*
*  Todd K. Moon
*  Utah State University
*
*  Date: 8/10/18
***************************************/

#ifndef POLARSTACK_H
#define POLARSTACK_H

#include <iostream>
using namespace std;


class polarstack {
public:
   int *stackmem;
   int stackptr;
   int maxstacksize;
   polarstack(void) {
	  stackmem = 0; // default: no space allocated yet
	  int stackptr = 0;
	  maxstacksize = 0;
   }
   polarstack(int L) {
	  stackmem = 0;
	  initspace(L);
   }
   ~polarstack() {
	  if(stackmem) delete[] stackmem;
   }
   void initspace(int L) {
	  if(stackmem == 0) {
	  	 stackmem = new int[L];
	  }
	  else {
	  	 delete [] stackmem;
	  	 stackmem = new int[L];
	  }
	  maxstacksize = L;
	  stackptr = 0;
   }

   void push(int val) {
	  if(stackptr == maxstacksize) {
		 std::cerr << "Error: polarstack full on push\n";
		 exit(-1);
	  }
	  stackmem[stackptr++] = val;
   }

   int pop() {
	  if(stackptr == 0) {
		 std::cerr << "Error: popping from empty stack\n";
		 exit(-1);
	  }
	  return stackmem[--stackptr];
   }
   void clearstack() {
	  stackptr = 0;
   }
   bool empty() {
	  if(stackptr == 0) return true;
	  else return false;
   }
   void printstack() {  // for testing/display purposes only
	  //	  cout << "number in stack=" << stackptr << ": ";
	  for(int i = 0; i < stackptr; i++) {
		 cout << stackmem[i] << " ";
	  }
	  cout << endl;
   }
};




#endif

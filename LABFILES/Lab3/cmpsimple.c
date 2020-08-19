/**************************************
*
*  Program: cmpsimple --- a simple (and slow) file comparison program
*
*
*  Todd K. Moon
*  Date: August 19. 2004
*
***************************************/

/* Copyright 2004 by Todd K. Moon
 Permission is granted to use this program/data
 for educational/research only
*/

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
   FILE *fout1,*fout2;
   unsigned char b1, b2;
   unsigned int byteno;

   if(argc==1) {
	  printf("Usage: %s file1 file2\n",argv[0]);
	  exit(-1);
   }
   fout1 = fopen(argv[1],"rb");
   fout2 = fopen(argv[2],"rb");
   byteno = 0;
   while(fread(&b1,1,1,fout1)) {
	  byteno++;
	  if(fread(&b2,1,1,fout2)) {
		 if(b1 != b2) {
			printf("Different bytes at byte number %ud\n",byteno);
		 }
	  }
	  else {
		 printf("Different file length\n");
		 break;
	  }
   }
   if(fread(&b2,1,1,fout2)) {
	  printf("Different file lengths\n");
   }
   fclose(fout1);
   fclose(fout2);
}
   

/*
Local Variables:
compile-command:"gcc -o cmpsimple cmpsimple.c"
End:
*/

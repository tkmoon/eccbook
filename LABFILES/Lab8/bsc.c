/**************************************
*
*  Program: bsc --- pass a file of data 
*  through a "binary symmetric channel"
*
*
*  Todd K. Moon
*  Date: Dec. 19, 2002
*
***************************************/

/* Copyright 2004 by Todd K. Moon
 Permission is granted to use this program/data
 for educational/research only
*/

#include <stdio.h>
#include <stdlib.h>
#define BUFLEN 100

int main(int argc, char *argv[])
{
   double p,p1;
   char *infname;
   char *outfname;
   FILE *fin, *fout;
   int i,j,nread,nwrote;
   unsigned char b;
   unsigned char databuf[BUFLEN];
   int K=0;						/* set if blockreport argument used */
   int bytecount = 0, blockcount = 0;
   int byteerr;
   int biterrcount = 0, byterrcount=0;;

   if(argc==1) {
	  printf("Usage: %s p infile outfile [blockreport]\n",argv[0]);
	  exit(-1);
   }
   p = atof(argv[1]);
   infname = argv[2];
   outfname = argv[3];
   if(argc==5)  {
	  K = atoi(argv[4]);
   }
   fin = fopen(infname,"r");
   fout = fopen(outfname,"wb");

   while((nread = fread(databuf,1,BUFLEN,fin))) { /* read a bunch of bytes */
	  for(i = 0; i < nread; i++) { /* work through each byte */
		 byteerr=0;
		 for(j = 0; j < 8; j++) { /* and each bit in the byte */
			p1 = (double)rand()/(double)RAND_MAX;  /* uniform [0,1] random */
			if(p1<p) {  /* bit error occurred */
			   databuf[i] ^= (1<<j); /* flip the bit */
			   biterrcount++;
			   byteerr = 1;
			}
		 }
		 if(byteerr) byterrcount++;
		 if(K) {
			bytecount++;
			if(bytecount % K == 0) {  /* a new block is starting */
			   printf("Block: %d   Bit Errors: %d  Byte Errors: %d\n",
					  blockcount,biterrcount,byterrcount);
			   biterrcount = byterrcount = 0;
			   blockcount++;
			}
		 }

	  }
	  nwrote = fwrite(databuf,1,nread,fout);
	  if(nwrote != nread) {
		 printf("Error: Cannot write to output file\n");
		 exit(-1);
	  }
   }
   fclose(fin);
   fclose(fout);
}
   

/*
Local Variables:
compile-command:"gcc -o bsc bsc.c"
End:
*/

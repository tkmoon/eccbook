int gcd(int a, int b)
/* A simple example of the Euclidean algorithm: */
/* Compute g = (a,b), where a>0 and b>0 */
/* (A better function would include sign checking) */

/* Copyright 2004 by Todd K. Moon
 Permission is granted to use this program/data
 for educational/research only
*/

{  int g;
   while(b) {
      g = b;
      b = a % b;  /* compute remainder of a/b */
      a = g;
   }
   return a;
}

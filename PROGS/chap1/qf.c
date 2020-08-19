/**********************
*
*  Computes Q(x) = \int_x^\infty \exp(-t*t/2)/\sqrt(2*\pi)
*  using a rational approximation
*  (See Abramowitz  Stegun, Handbook of Mathematical Functions,
*    p 932, formula 26.2.17)
*
**************************/
/* Copyright 2004 by Todd K. Moon
 Permission is granted to use this program/data
 for educational/research only
*/

#include <math.h>

double qf(double x)
{
double t,dq;
double q;

   if(x < -7) return(1.0);
   t=1.0/(1.0+0.2316419*fabs(x));
   dq = ((((1.330274429*t - 1.821255978)*t + 1.78147793)*t - .356563782)*t
   + .319381530)*t * exp(-x*x*0.5) * .39894228;
   if(x > 0.0)
      q = dq;
   else
      q = (1.0-dq);
   return(q);
}

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  
#include <nicklib.h>  

#include "kimsubs.h" 

static int verbose = NO ;
static int debug = 0 ;

#define N 1000 

/** 
 simply prinbt Kimura function for given x, t 
*/
int main  (int argc , char **argv) 
{
   double x, t, y, ktrans ;
   int j ;
   
  x = atof(argv[1]) ;
  t = atof(argv[2]) ;
  printf("## x:   %10.4f  t:  %10.4f\n", x, t) ;
  for (j=0; j<=N; j++) { 
   y = (double) j ;  y /= (double) N ;
   ktrans = evalkim(y, x, t, 1.0) ;
   printf("%9.3f  %12.6f\n", y, ktrans) ;
  }
  return 0 ;
}

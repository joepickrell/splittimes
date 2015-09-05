#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>
#include <getpars.h> 


static int verbose = NO ; 
#define  VERSION  "100"  

int modehprob(int n, int a, int m)  ;

int main() 
{
  int n=50,  m=16, a=5  ;
  int k, t, iter ;
  int hist[100] ;
  double y, y1, y2  ;

  ivzero(hist, 100) ;
  for (iter =1 ; iter <= 10000; ++iter) { 
   n = 63 ;  
   m = ranmod(n-1) + 1 ; 
   a = ranmod(n-1) + 1 ; 
   
   t = modehprob(n, a, m) ;  
   y1 = loghprob(n, a, m, t) ;
   y2 = loghprob(n, a, m, t+1) ;
   if (y2>y1) printf("hit: %d %d %d %d %12.6f %12.6f\n", n, a, m, t, y1, y2) ;
   y2 = loghprob(n, a, m, t-1) ;
   if (y2>y1) printf("hitx: %d %d %d %d %12.6f %12.6f\n", n, a, m, t, y1, y2) ;
  }

}


int modehprob(int n, int a, int m) 
{
 double v, y ;
 double pm, logpm, w, ru, rw, rat ;
 int mode, k, x, zans ;

 v = (double) (a+1)*(m+1) / (double) (n+1) ; 
 mode = (int) v ;
 for (;;) { 
  if (loghprob(n, a, m, mode) < loghprob(n, a, m, mode+1)) { 
   mode = mode + 1 ;
   continue ;
  }
  if (loghprob(n, a, m, mode) < loghprob(n, a, m, mode-1)) { 
   mode = mode - 1 ;
   continue ;
  }
  break ; 
 }

 return mode ; 

}

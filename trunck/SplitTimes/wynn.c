#include <nicklib.h>  
#include <math.h>  

#define N 20 

double wynn(double *v, int n, double *acc, int *nacc)   
{
 double *x0, *x1, *xn ; 
 double t, amax, amin ;
 int iter = 0, j, nn  ;

 vmaxmin(v, n, &amax, &amin) ;  
 if (amax<=amin) {  
  vclear(acc, amax, n/2) ;
  *nacc = n/2 ;
  return amax ;
 }

 ZALLOC(x0, n, double) ; 
 ZALLOC(x1, n, double) ;
 ZALLOC(xn, n, double) ;
 copyarr(v, x1, n) ;  
 nn = n ;  
 for (;;) {  
  for (j=0; (j+1) < nn ; ++j) {  
   t = x0[j+1] + 1.0/(x1[j+1]-x1[j]) ;
   xn[j] = t ;
  }
  --nn ; 
  if (nn<2) break ;  
  copyarr(x1, x0, n) ;
  copyarr(xn, x1, n) ;

  for (j=0; (j+1) < nn ; ++j) {  
   t = x0[j+1] + 1.0/(x1[j+1]-x1[j]) ;
   xn[j] = t ;
  }
  --nn ; 
  if (nn<2) break ;  
  copyarr(x1, x0, n) ;
  copyarr(xn, x1, n) ;
  acc[iter] = t ; 
  ++iter ;
 }
 free(x0) ; free(x1) ; free(xn) ;
 *nacc = iter ;
 return t ;
}


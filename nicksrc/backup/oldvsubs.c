#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "strsubs.h" 
#include "vsubs.h" 

/** 
 tiny routines BLAS? 
 a small library to do simple arithmetic
 on 1D vectors with no skips 
*/
void 
vsp(double *a, double *b, double c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] + c ;
}
void 
vst(double *a, double *b, double c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] * c ;
}
void 
vvt(double *a, double *b, double *c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] * c[i] ;
}
void 
vvp(double *a, double *b, double *c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] + c[i] ;
}
void 
vvm(double *a, double *b, double *c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] - c[i] ;
}
void 
vvd(double *a, double *b, double *c, int n) 
{
   int i ;
   for (i=0; i<n; i++)  {
    if (c[i] == 0.0) 
      fatalx("(vvd): zero value in denominator\n") ;
    a[i] = b[i] / c[i] ;
   }
}
void 
vsqrt(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    if (b[i]<0.0) 
     fatalx("(vsqrt): negative value %g\n",b[i]) ;
    if (b[i] == 0.0) {
     a[i] = 0.0 ;
     continue ;
    }
    a[i] = sqrt(b[i]) ;
   }
}
void 
vinvert(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    if (b[i] == 0.0) 
     fatalx("(vinvert): zero value\n") ;
    a[i] = 1.0 / b[i] ;
   }
}
void 
vabs(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    a[i] = fabs(b[i]) ;
   }
}
void 
vlog(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    if (b[i]<=0.0) 
     fatalx("(vlog): negative or zero value %g\n",b[i]) ;
    a[i] = log(b[i]) ;
   }
}
void 
vlog2(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    if (b[i]<=0.0) 
     fatalx("(vlog2): negative or zero value %g\n",b[i]) ;
    a[i] = log2(b[i]) ;
   }
}

void 
vexp(double *a, double *b,  int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    a[i] = exp(b[i]) ;
   }
}
void 
vclear(double *a,  double c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] =  c ;
}
void 
vzero(double *a,  int n) 
{
  vclear(a, 0.0, n) ;
}
void 
cclear(char *a,  char c, int n) 
/** 
 be careful nothing done about NULL at end
*/
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] =  c ;
}
void 
ivvp(int *a, int *b, int *c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] + c[i] ;
}
void 
ivsp(int *a, int *b, int c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] = b[i] + c ;
}
void 
ivclear(int *a,  int c, int n) 
{
   int i ;
   for (i=0; i<n; i++) 
    a[i] =  c ;
}
void
ivzero(int *a,  int n) 
{
  ivclear(a, 0, n) ;
}
void 
vclip(double *a, double *b,double loval, double hival,int n)  
{
/* clip off values to range [loval,hival] */
   int i ;
   double t ;

   for (i=0; i<n; i++) {
    t = MAX(b[i],loval) ;
    a[i] = MIN(t,hival) ;
   }
}
void vmaxmin(double *a, int n, double *max, double *min)  
{

    int i ;
    double tmax, tmin ;

    tmax = tmin = a[0] ;
    for (i=1; i<n; i++) {
      tmax = MAX(tmax, a[i]) ;
      tmin = MIN(tmin, a[i]) ;
    }
    if (max != NULL) *max = tmax ;
    if (min != NULL) *min = tmin ;
}
void vlmaxmin(double *a, int n, int *pmax, int *pmin)  
/** 
 return location 
*/
{

    int i ;
    double tmax, tmin ;
    double lmax, lmin ;

    tmax = tmin = a[0] ;
    lmax = lmin = 0 ;
    for (i=1; i<n; i++) {
      if (a[i]>tmax) {
       tmax = a[i] ;
       lmax=i ;
      }
      if (a[i]<tmin) {
       tmin = a[i] ;
       lmin=i ;
      }
    }
    if (pmax != NULL) *pmax = lmax ;
    if (pmin != NULL) *pmin = lmin ;
}
void ivmaxmin(int *a, int n, int *max, int *min)  
{

    int i ;
    int tmax, tmin ;

    tmax = tmin = a[0] ;
    for (i=1; i<n; i++) {
      tmax = MAX(tmax, a[i]) ;
      tmin = MIN(tmin, a[i]) ;
    }
    if (max != NULL) *max = tmax ;
    if (min != NULL) *min = tmin ;
}
void ivlmaxmin(int *a, int n, int *pmax, int *pmin)  
/** 
 return location 
*/
{

    int i ;
    int tmax, tmin ;
    int lmax, lmin ;

    tmax = tmin = a[0] ;
    lmax = lmin = 0 ;
    for (i=1; i<n; i++) {
      if (a[i]>tmax) {
       tmax = a[i] ;
       lmax=i ;
      }
      if (a[i]<tmin) {
       tmin = a[i] ;
       lmin=i ;
      }
    }
    if (pmax != NULL) *pmax = lmax ;
    if (pmin != NULL) *pmin = lmin ;
}
double 
vdot(double *a, double *b, int n) 
{
   int i; 
   double ans=0.0 ;
   for (i=0; i<n; i++) 
    ans += a[i]*b[i] ;

   return ans ;

}
void 
getdiag(double *a, double *b, int n)  
/* extract diagonal */
{
   int i, k ;

   for (i=0; i<n; i++) {
    k = i*n+i ;
    a[i] = b[k] ;
   }
}

void copyarr(double *a,double *b,int n) 
{
  int i ;
  for (i=0; i<n; i++) {
   b[i] = a[i] ;
  }
}

void copyiarr(int *a,int *b,int n) 
{
  int i ;
  for (i=0; i<n; i++) {
   b[i] = a[i] ;
  }
}
void copyiparr(int **a,int **b,int n) 
{
  int i ;
  for (i=0; i<n; i++) {
   b[i] = a[i] ;
  }
}
void dpermute(double *a, int *ind, int len) 
{

  int i , k ;
  double *rrr ;

  ZALLOC(rrr, len, double) ;

  for (i=0; i<len; i++) {
   rrr[i] = a[i] ;
  }

  for (i=0; i<len; i++) {
   k = ind[i] ;
   a[i] = rrr[k] ;
  }

  free (rrr) ;
}
void ipermute(int *a, int *ind, int len) 
{

  int i , k ;
  int *rrr ;

  ZALLOC(rrr, len, int) ;

  for (i=0; i<len; i++) {
   rrr[i] = a[i] ;
  }

  for (i=0; i<len; i++) {
   k = ind[i] ;
   a[i] = rrr[k] ;
  }

  free (rrr) ;
}
void dppermute(double **a, int *ind, int len) 
{

  int i , k ;
  double **rrr ;

  ZALLOC(rrr, len, double *) ;

  for (i=0; i<len; i++) {
   rrr[i] = a[i] ;
  }

  for (i=0; i<len; i++) {
   k = ind[i] ;
   a[i] = rrr[k] ;
  }

  free (rrr) ;
}
void ippermute(int **a, int *ind, int len) 
{

  int i , k ;
  int **rrr ;

  ZALLOC(rrr, len, int *) ;

  for (i=0; i<len; i++) {
   rrr[i] = a[i] ;
  }

  for (i=0; i<len; i++) {
   k = ind[i] ;
   a[i] = rrr[k] ;
  }

  free (rrr) ;
}

double  asum(double *a, int n) 
{
   int i; 
   double ans=0.0 ;
   for (i=0; i<n; i++) 
    ans += a[i] ;

   return ans ;
}
int  intsum(int *a, int n) 
{
   int i; 
   int ans=0 ;
   for (i=0; i<n; i++) 
    ans += a[i] ;

   return ans ;
}


double  aprod(double *a, int n) 
/* overflow not checked */
{
   int i; 
   double ans=1.0 ;
   for (i=0; i<n; i++) 
    ans *= a[i] ;

   return ans ;
}
double  asum2(double *a, int n) 
{
   int i; 
   double ans=0.0 ;
   for (i=0; i<n; i++) 
    ans += a[i]*a[i] ;

   return ans ;
}

double trace(double *a, int n) 
{
   double *diags, t ;
   ZALLOC(diags,n,double) ;
   getdiag(diags,a,n) ; /* extract diagonal */
   t = asum(diags,n) ;
   free(diags) ;
   return t ;
}

int nnint(double x)
{
   return (int) rint(x) ;
}
void 
countcat(int *tags, int n,int *ncat,int nclass)  
/* simple frequency count of integer array */

{
   int i, k; 
   ivclear(ncat, 0, nclass) ;
   for (i=0 ; i<n ; i++)  {
    k = tags[i] ;
    if ( (k<0) || (k >= nclass)) 
     fatalx("(countcat) bounds error\n") ;
    ++ncat[k] ;
   }
}
void rowsum(double *a, double *rr, int n) 
{
   int i,j ;
   vclear(rr,0.0,n) ;
   for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
     rr[j] += a[i+j*n] ;
    }
   }
}
void colsum(double *a, double *cc, int n) 
{
   int i,j ;
   vclear(cc,0.0,n) ;
   for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
     cc[i] += a[i+j*n] ;
    }
   }
}
void rrsum(double *a, double *rr, int m, int n) 
{
   int i,j ;
   vclear(rr,0.0,n) ;
   for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
     rr[j] += a[i+j*m] ;
    }
   }
}
void ccsum(double *a, double *cc, int m, int n) 
{
   int i,j ;
   vclear(cc,0.0,m) ;
   for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
     cc[i] += a[i+j*m] ;
    }
   }
}
void printmat(double *a, int m, int n) 
/** 
 print a matrix n wide m rows  
*/
{
      int i,j, jmod ;
      for (i=0; i<m; i++)  {
        for (j=0; j<n; j++)  {
         printf("%9.3f ", a[i*n+j]) ;
         jmod = (j+1) %5 ;
         if (jmod == 0) {
          printf("  ...\n") ;
         }
        }
        printf("\n") ;
      }
}
void floatit(double *a, int *b, int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    a[i] = (double) b[i] ;
   } 
}

void fixit(int  *a, double *b, int n) 
{
   int i ;
   for (i=0; i<n; i++) {
    a[i] = nnint(b[i]) ;
   } 
}
int findfirst(int *a, int n, int val) 
{
   int i ;
   for (i=0 ; i<n; i++)  {
    if (a[i] == val) return i ;
   }
   return -1 ;  
}
void idperm(int *a, int n) 
{
     int i ;
     for (i=0; i<n; i++) 
      a[i] = i ;
}
double log2(double y) 
{
  if (y<=0.0) fatalx("(log2) negative argument\n") ;
  return (log(y)/log(2.0)) ;
}
double log2fac(int n) 
/** 
 log base2 (factorial n))
*/
{
  double y, x ;
  x = (double) (n+1) ;
  y = lgamma(x) ;
  return (y/log(2.0)) ;
}

double addlog(double a, double b) 
{
 /* given a = log(A)
          b = log(B)
    returns log(A+B) 
    with precautions for overflow etc
*/
    double x, y, z ;

    x = MIN(a,b) ;
    y = MAX(a,b) ;

/** 
 answer is log(1+A/B)+log (B)  
*/
    z = x-y ;  
    if (z < -50.0)  return y ;
    z = 1.0+exp(z) ;
    z = log(z) + y ;
    return (z) ;

}



double vldot(double *x, double *y, int n) 
/** 
 x. log(y) 
*/
{
  double *z, ans ;
  ZALLOC(z, n, double) ;
  vsp(z, y, 1.0e-20, n) ;
  vlog(z, z, n) ;
  ans = vdot(x, z, n) ;
  free (z) ;
  return ans ;
}

double pow10 (double x) 
{
     return exp(x*log(10.0)) ;
}


double vpow10 (double *a, double *b, int n) 
{
     int i ;
     for (i=0; i< n; i++)  
      a[i] = exp(b[i] * log(10.0)) ;
}

double vlog10 (double *a, double *b, int n) 
{
     int i ;
     for (i=0; i< n; i++)  
      a[i] = log10(b[i]) ;
}


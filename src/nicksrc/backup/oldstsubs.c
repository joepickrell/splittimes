#include <stdio.h>
#include <math.h> 

#include "statsubs.h" 
#include "vsubs.h" 

#define EPS1 .001 
#define EPS2 1.0e-8 
#define ZLIM 12 
#define QSIZE 100 

static double zzprob(double zval) ;
static int firstgt(double val, double *tab, int n) ;

double ntail(double zval) 
/* uses erfc */
{
  double pi,  t ;
  double p, q, d ;

  if (zval<ZLIM) {
  t = zval/sqrt(2.0) ;
  q = erfc(t)/2.0 ;
  return q ;
  }

  pi = 2.0*acos(0.0)  ;

  t = exp(-0.5*zval*zval) ;
  t /= (sqrt(2.0*pi) * zval) ;

  return t ;

}

double zzprob(double pval) {
  double x, dev, p, q, d, h, u ;
  double pi ;
   int iter ;

/** 
 approximate normal by 1/(sqrt 2 pi) * exp (-0.5*x*x) / x  
 appropriate for pval very small
*/  
/* Feller I  page 166 */

  if (pval==0.0) return 50.0 ;

  pi = 2.0*acos(0.0) ;
  u = -log(sqrt(2.0*pi)*pval) ;
/* solve x*x/2 + log(x) = u */

  x = sqrt(2.0*u) ;
  for (iter=1; iter<=10; ++iter) {  
   q = (0.5*x*x) + log(x) ;
   d = x + (1.0/x) ;
   dev = u - q;
   h = dev/d ;
/**
   printf("zz %12.6e  %12.6e  %12.6e\n",x,dev,h) ;
*/
   x += h ;
   if (fabs(h)<1.0e-7) return x ;
  }
  return x ;
}

double medchi(int *cls, int len, int *n0, int *n1, double *xtail) 
{
/* compute 2x2 chisq splitting at median */
 int i, m0,m1,n,m ;
 double arr[4], y, ys, p, q, d ;
 *n0 = *n1 = 0 ;
 for (i=0; i<len; i++) {
  if (cls[i]>1) continue ;
  if (cls[i]<0) continue ;
  if (cls[i]==0) ++*n0 ;
  if (cls[i]==1) ++*n1 ;
 }
 if (MIN(*n0,*n1)==0) {
/**
  for (i=0; i<len ; i++) {
   printf("zz1 %d %d\n",i,cls[i]) ;
  }
*/
  *xtail = 1.0 ;
  return 0 ;
 }
 m = (*n0+*n1)/2 ;
 m0 = m1 = 0;
 for (i=0; i<len; i++) {
  if (cls[i]>1) continue ;
  if (cls[i]<0) continue ;
  if (cls[i]==0) ++m0 ;
  if (cls[i]==1) ++m1 ;
  if ((m0+m1) == m)  break ;
 }

 arr[0] = (double) m0 ;
 arr[1] = (double) m1 ;
 arr[2] = (double) (*n0-m0) ;
 arr[3] = (double) (*n1-m1) ;

 y = conchi(arr,2,2) ;
 ys = sqrt(y+EPS2) ;
 q = ntail(ys) ;
 
 *xtail = q ;

/**
 printf("zzchsq ") ;
 for (i=0; i<=3 ; i++) 
  printf("%9.3f ",arr[i]) ;
 printf("%9.3f %9.3f\n",y,q) ;
*/

 return y ;

}

double ks2(int *cls, int len, int *n0, int *n1, double *kstail)
{
/*
 compute KS statistic 
 cls should be 0 or 1.  if larger take as invalid
*/
 int i ;
 double en0, en1, en ; 
 double y, ymax  ;
/* count class sizes */

 if (len <= 1) {
  *kstail = 1.0 ;
  return 0 ;
 }

 *n0 = *n1 = 0 ;
 for (i=0; i<len; i++) {
  if (cls[i]>1) continue ;
  if (cls[i]<0) continue ;
  if (cls[i]==0) ++*n0 ;
  if (cls[i]==1) ++*n1 ;
 }
 if (MIN(*n0,*n1)==0) {
/**
  printf("warning ks2 has only 1 class passed\n") ;
  for (i=0; i<len ; i++) {
   printf("zz1 %d %d\n",i,cls[i]) ;
  }
*/
  *kstail = 1.0 ;
  return 0 ;
 }

 en0 = (double) *n0 ;
 en1 = (double) *n1 ;
 

 ymax = y = 0.0 ;  /* running stat */ ;
 for (i=0; i<len; i++) {
  if (cls[i]>1) continue ;
  if (cls[i]<0) continue ;
  if (cls[i]==0) y += 1.0/en0 ;
  if (cls[i]==1) y -= 1.0/en1 ;
  ymax = MAX(ymax,fabs(y)) ;
 }

/*  Numerical recipes p 626 */
 en = sqrt(en0*en1/(en0+en1)) ;
 y = en+.12+(0.11/en) ;
 y *= ymax ;
 *kstail = probks(y) ;
 return y ;
/** crude analysis:  
  variance of 1 step above is (1/n0 + 1/n1) / (n0+n1) 
  and so variance of y is brownian motion not bridge is (1/n0+1/n1) 
  We want to rescale y to correspond to Brownian bridge.  
  First order correction is en.  We actually use 
   a Bartlett correction of some sort 
  Normalized y seems like what to return.
*/

}
double probks(double lam) 
/* KS tail area: Numerical recipes p 626 */
{
  int j ;
  double a2, fac=2.0, sum=0.0, term, termbf=0.0 ;
  double t ;

  a2 = -2.0*lam*lam ;
  for (j=1; j<=100; j++) {
   t  = a2* (double) (j*j) ;
   term = fac*exp(t) ;
   sum += term ;
   t = fabs(term) ;
   if ((t <= EPS1*termbf) || ( t <= EPS2*sum)) return sum ;
   fac = -fac ;
   termbf = fabs(term) ;
  }
  return 1.0 ;
}
double conchi(double *a, int m, int n) 
/* a is m rows n columns.  contingency chisq */
{
 double *rsum, *csum, ee, tot=0, chsq=0, y ;
 int i,j,k ;

 ZALLOC(rsum,m,double) ;
 ZALLOC(csum,n,double) ;

 for (i=0; i<m; i++) {  
  for (j=0; j<n; j++) {  
   k = i*n+j ;
   rsum[i] += a[k] ;
   csum[j] += a[k] ;
   tot += a[k] ;
  }
 }
 if (tot < 0.001) 
   fatalx("(conchi) no data\n") ;
 for (i=0; i<m; i++) {  
  for (j=0; j<n; j++) {  
   k = i*n+j ;
   ee = rsum[i]*csum[j]/tot ;
   if (ee < EPS2) 
    fatalx("(conchi) zero row or column sum\n") ;
   y = a[k]-ee ;
   chsq += (y*y)/ee ;
  }
 }
 free(rsum) ;  free(csum) ;
 return chsq ;
}

double chitest(double *a, double *p, int n) 
/* a is n boxes.  Goodness of fit test to p */
{
 
 double *x, *b, *pp ;
 double y1=0.0, y2=0.0 ;
 int i ;

 ZALLOC(pp, n, double) ;
 if (p != NULL)
  copyarr(p,pp,n) ;
 else 
  vclear(pp, 1.0, n) ;

 y1 = asum(pp,n) ;
 y2 = asum(a,n) ;

 if ( (y1==0.0) || (y2==0.0) ) { 
  free(pp) ;
  return 0.0 ;
 }

 ZALLOC(x,n,double) ;
 ZALLOC(b,n,double) ;


 vst (x, pp, y2/y1, n) ;  /* expected */

 vsp (x, x, .0001, n) ;
 vvm (b, a, x, n) ;  
 vvt (b, b, b, n) ;
 vvd (b, b, x, n) ;

 y1 = asum(b,n) ;

 free(x) ;
 free(b) ;

 return y1 ;

}
double zprob(double ptail) 
/** inverse normal */
/** 
 This crude routine could be easily modified 
 to return result to machine accuracy, at least for 
 z < ZLIM (use a Newton iteration) but for likely 
 uses this is irrelevant.

 Method: 
1) Build table 0..ZLIM (1/QSIZE) of normal tail.
2) Do linear interpolation
*/
 
 
{
#define QSIZE 100
  static int first = YES ;
  static double ptiny ;
  static int numbox = QSIZE*ZLIM ;
  static double *ztable, *ptable ;
  double z, p, t, plog; 
  double ylo, yhi, ya, yb ; 
  int i, k ;
  if (first=YES) {  
   ZALLOC(ptable, numbox, double) ;
   ZALLOC(ztable, numbox, double) ;
   first = NO ;
   for (z=0.0, i=0; i<numbox; i++) {
    p = ntail(z) ;
    ztable[i] = z ;  
    ptable[i] = -log(p) ;
    z += 1.0/(double) QSIZE ;
  }
  ptiny = p ;
 }
 if (ptail==0.5) return 0.0 ;
 if (ptail>0.5) {
  z = zprob(1.0-ptail) ;
  return -z ;
 }
 if (ptail<ptiny) {
  z = zzprob(ptail) ;
  return z ;
 }
/** replace by binary or interpolating search */  
 plog = -log(ptail) ;
 k = firstgt(plog, ptable, numbox) ;
 if (k==0) return ztable[0] ;
 if (k==numbox) return ztable[numbox-1] ;
 ylo = ptable[k-1] ;
 yhi = ptable[k] ;
 ya = (yhi-plog) ;
 yb = plog-ylo ;
 z = (ya*ztable[k-1]+yb*ztable[k])/(ya+yb) ;
 return z ;
}
int firstgt(double val, double *tab, int n) 
{
/* tab sorted in ascending order */
  int i ;

  if (val>tab[n-1]) return n ;
  for (i=0; i<n; i++) {
   if (val<=tab[i]) return i ;
  }
}


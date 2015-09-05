#include <stdio.h>
#include <math.h> 

#include "statsubs.h" 
#include "vsubs.h" 

#define EPS1 .001 
#define EPS2 1.0e-12 
#define ZLIM 20 
#define QSIZE 10 


static double *bern ; /* bernouilli numbers */
static int bernmax = 0 ;

static double zzprob(double zval) ;
static double znewt(double z, double ptail) ;

double nordis(double zval) 
/* normal density */
{
  double pi,  t ;

  pi = 2.0*acos(0.0)  ;

  t = exp(-0.5*zval*zval) ;
  t /= sqrt(2.0*pi)  ;

  return t ;

}

double ntail(double zval) 
/** normal distribution tail area 
 uses erfc 
*/

{
  double pi,  t ;
  double p, q, d ;

  if (zval == 0.0) return 0.5 ;
  if (zval<0.0) return (1.0 - ntail(-zval)) ;
 
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

/** approximate normal by 1/(sqrt 2 pi) * exp (-0.5*x*x) / x   */  
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
   if (ee < EPS2) {
    printf("bad conchi\n") ;
    printmat(a, m, n) ;
    fatalx("(conchi) zero row or column sum\n") ;
   }
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
{
  static int first = YES ;
  static double ptiny ;
  static int numbox = QSIZE*ZLIM ;
  static double *ztable, *ptable ;
  double z, p, t, plog; 
  double ylo, yhi, ya, yb ; 
  int i, k ;
  if (first==YES) {  
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
  return znewt(z,ptail) ;
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
 return znewt(z, ptail) ;
}
double znewt(double z, double ptail) 
{
/** 
 newton step 
 z is very good approximation 
*/
    double p0, pder, h ;
    double pi, zz ;
    int iter ;
    pi = 2.0*acos(0.0) ;
    zz = z ;
    for (iter = 1; iter<=2; ++iter) {
     p0 = ntail(zz) ;
     pder = -exp(-0.5*zz*zz)/sqrt(2*pi) ;
     if (pder==0.0) return zz ;
     h = (ptail-p0)/pder ;
     if (h==0.0) return zz ;
     zz += h ;
    }
    return zz ;
}
int ifirstgt(int val, int *tab, int n) 
{
/* tab sorted in ascending order */
  int i ;

  if (val>=tab[n-1]) return n ;
  for (i=0; i<n; i++) {
   if (val<tab[i]) return i ;
  }
}

int firstgt(double val, double *tab, int n) 
{
/* tab sorted in ascending order */
  int i ;

  if (val>=tab[n-1]) return n ;
  for (i=0; i<n; i++) {
   if (val<tab[i]) return i ;
  }
}

void mleg(double a1, double a2, double *p, double *lam) 
{
   int iter ;
   double s, pp, ll  ;
   double top, bot, fval ;
   int debug = NO ;

/** 
 solve 
 p/lam = a1 ; psi(p) - log(lam) = a2 ;
 Thus psi(p) - log(p) = a2 - log(a1) 
*/
  s = a2 - log(a1) ;   

  if (s>=0.0) fatalx("log E(x) < E(log (x)) \n") ;
  pp = -s ;

  for (iter = 1; iter <= 30 ; ++iter) {  
   fval = s - (psi(pp) - log (pp)) ;
   if (debug)
    printf("yy1 %3d %9.3f %9.3f\n",iter,pp,fval) ;
   if (fval<0.0)  break ;
   pp *= 2.0 ;
  }

  for (iter = 1; iter <= 30 ; ++iter) {  
   fval = s - (psi(pp) - log (pp)) ;
   if (fval>0.0)  break ;
   if (debug)
    printf("yy2 %3d %9.3f %9.3f\n",iter,pp,fval) ;
   pp /= 2.0 ;
  }
  
  for (iter = 1; iter <= 10 ; ++iter) {  
   fval = psi(pp) - log (pp) ;
   top = s-fval ;
   bot =  tau(pp) - (1.0/pp) ;
   if (debug)
    printf("%3d %12.6f %12.6f\n",iter,pp,top) ;
   pp += top/bot ;
  }
   ll = pp/a1 ;
   *p = pp  ;
   *lam = ll ;
}


double psi(double x) 
{
 double y, zz, term ;
 int k ;

 bernload() ; 
 if (x<10.0) return (psi(x+1.0) - 1.0/x) ;
 y = log(x) - 1.0/(2.0*x) ;
 zz = 1.0 ;
 for (k=1; k<= bernmax/2 ; k++)  {
  zz /= (x*x) ;
  term = bernum(2*k)/(double) (2*k) ;
  term *= zz ;
  y -= term ;
 }
 return y ;
}
double tau(double x) 
/*
 derivative of psi 
*/
{
 double y, zz, term ;
 int k ;

 bernload() ;
 if (x<10.0) return (tau(x+1.0) + 1.0/(x*x)) ;

 y = 1.0/x  + 1.0/(2.0*x*x) ;
 zz = 1.0/x ;
 for (k=1; k<= bernmax/2 ; k++)  {
  zz /= (x*x) ;
  term = bernum(2*k)/(double) (2*k) ;
  term *= zz ;
  term *= - (double) (2*k) ;
  y -= term ;
 }
 return y ;
}
void bernload() 
{
 if (bernmax>0) return ;
 bernmax = 14 ;
 ZALLOC(bern, bernmax+1, double) ;
 bern[0] = 1.0 ;
 bern[1] = -1.0/2.0 ;
 bern[2] =  1.0/6.0 ;
 bern[4] =  -1.0/30.0 ;
 bern[6] =   1.0/42.0 ;
 bern[8] =  -1.0/30.0 ;
 bern[10] =  5.0/66.0  ;
 bern[12] = -691.0/2730.0 ;
 bern[14] =  7.0/6.0 ;
}
double bernum(int k) 
{
  bernload() ;
  if ((k<0) || (k>bernmax)) fatalx("(bernum) bad arg: %s\n",k) ;
  return (bern[k]) ;
}

double dilog(double x) 
{
 return li2(x) ;
}

double li2(double x) 
{
 double pi, sum=0.0, term, top, z ;
 int k ;

 pi = acos(0.0)*2.0 ;
 if (x<=0.0) return pi*pi/6.0 ;
 if (x>=1.0) return 0 ;
 if (x<0.5) {
  return (-log(x)*log(1-x) + (pi*pi/6.0) -li2(1.0-x)) ;
 }
  z = 1-x ;
  top = 1.0 ;
  for (k=1; k<= 100; k++) { 
   top *= z ;
   term = top/(double) (k*k) ;
   sum += term ;
   if (term <= 1.0e-20) break ;
  }
  return sum ;
}

double hwstat(double *x) 
/** Hardy-Weinberg equilibrium test 
    returns standard normal in null case. 
    +sign is excess heterozygosity 
    x[0] [1] [2] are counts for homozm hetero homo (alt allele)
*/
{

     double p, q, ysum, s1, y1, y2, ychi, sig ;
     double a1[3], a2[3] ;

     ysum = asum(x,3) ;
     if (ysum < 0.001) return 0.0 ;
     s1 = 2*x[2]+x[1] ;
     p  = 0.5*s1/ysum;
     q = 1.0-p ;

     a1 [0] = q*q ;
     a1 [1] = 2*p*q ;
     a1 [2] = p*p ;

     vsp(a1, a1, 1.0e-8, 3) ;
     vst(a2, x, 1.0/ysum, 3) ;
     vsp(a2, a2, 1.0e-8, 3) ;

     y2 = vldot(x, a2, 3) ;
     y1 = vldot(x, a1, 3) ;

     ychi = 2.0*(y2-y1) ;
     sig = sqrt(ychi+1.0e-8) ;

     if (a2[1]<a1[1]) sig = -sig ;
/* negative => hets lo */

     return sig ;
   

}

double bprob(double p, double a, double b) 
{
   double q, yl ;
   int iter ;

   q = 1.0 - p ;

   for (iter = 1; iter <= 10; ++iter) {
    yl = (a-1) * log(p) + (b-1) * log (q) ;
    yl -= lbeta(a, b) ;
    if (finite(yl)) return yl ;
// absurd hack.  NJP 
    if (p<=0.0) break ;
    if (q<=0.0) break ;
   }
   fatalx("bad bprob %12.3e %12.3e %12.3e\n", p, a, b) ;
}

double lbeta(double a, double b) 
{
   return (lgamma(a) + lgamma(b) - lgamma(a+b) ) ;
}


double dawson(double t) 
/** 
 Dawson's Integral 
 [A + S 7.31]
 exp(-t*t) \int ( exp(x^2), x = 0..t) 
 loosely based on mcerror.for 
*/
{
        
        double  z1,  cs, cr, cl ;
        double z1sq, cer ;
        int k ;

        z1 =  fabs(t) ;
        if (z1 <= 1.0e-8) return t ; 
/* derivative is 1 at 0 */
        z1sq = -t*t ;
        if (z1 < 4.5) {  
         cs = cr = z1 ;
         for (k= 1; k <= 999; ++k) {  
            cr *= z1sq/((double) k + 0.5) ;
            cs += cr ;
            if (fabs(cr/cs) < 1.0e-15) break ;
         }
         cer = cs ;
        }
        else {
         cl = 1/z1 ;
         cr = cl ;
         for (k=1; k<=13; ++k) {
          cr *= -((double) k-0.5) / z1sq ;
          cl += cr ;
          if (fabs(cr/cl) < 1.0e-15) break ;
         }
         cer = 0.5*cl ;
        }
        if (t<0) cer = -cer ; 
        return cer ;
}

double binlogtail(int n, int t, double p, char c) 
{

    double *bindis ;
    double val, base ;

    ZALLOC(bindis, n+1, double) ;
    genlogbin(bindis, n, p) ;
    base = bindis[t] ;
    vsp(bindis, bindis, -base, n+1) ;
    if (c=='+') {
     vexp(bindis+t, bindis+t, n-t+1) ;
     val = asum(bindis+t, n-t+1) ; 
    }
    else  {
     vexp(bindis, bindis, t) ;
     val = asum(bindis, t) ;
    }
    free(bindis) ;
    return (log(val) + base) ;

}

double binomtail(int n, int t, double p, char c) 
{
/** 
 c = '+':   P(S>=t) 
 c = '-':   P(S<t) 
 WARNING <= t use binomtail(n, t+1, ... 
*/
    double *bindis ;
    double val ;

    ZALLOC(bindis, n+1, double) ;
    genlogbin(bindis, n, p) ;
    vexp(bindis, bindis, n+1) ;
    if (c=='+') 
     val = asum(bindis+t, n-t+1) ; 
    else  
     val = asum(bindis, t) ;
    free(bindis) ;
    return val ;
}

void
genlogbin(double *a, int n, double p)
/* generate log prob for binomial distribution */
{
    double q, plog, qlog ;
    double *lfac, x ;
    int i, r, s ;
    q = 1.0-p ;
    plog = log(p), qlog = log(q) ;

    ZALLOC(lfac,n+1,double) ;
    for (i=1; i<=n; i++) {
     x = (double) i ;
     lfac[i] = lfac[i-1] + log(x) ;
    }
    for (r=0; r<=n; r++) {
     s = n-r ;
     x = lfac[n]-lfac[r]-lfac[s] ;  /* log binom coeff */
     x += ((double) r) * plog ;
     x += ((double) s) * qlog ;
     a[r] = x ; /* log prob */
    }
    free(lfac) ;
}


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  
#include <statsubs.h>  
#include <vsubs.h>  

#include "kimsubs.h" 
double wynn(double *v, int n, double *acc, int *nacc)   ;

// includes wynn acceleration

static int verbose = NO ;
static int debug = 0 ;

static int maxpol = MAXPOL ;

static double rfac3 (int m) ;
double intx2d(double y, int kdeg, double tt) ;

double evalkim(double x, double p, double t, double s2)  
{
  return  p * (1.0-p) * evalkim0(x, p, t, s2) ;
}

double evalkim0(double x, double p, double t, double s2) 
{

/** 
 evaluate forward transition prob for a 
 neutral Wright-Fisher diffusion 
 divided by p * (1-p) Makes sense as limit p->0 / p 
*/

   double tt  ;
   double xsum[MAXPOL], acc[MAXPOL] ;
   double sum = 0.0, term, alpha, q, z, xpol  ;
   int j, nacc ;
   static int first = YES ;

   if (p != pco) {  
    if (x == pco) return evalkim0(p,x,t,s2) ;
    setccoeff(p) ;
   }
   vzero(xsum, maxpol) ;
   z = 1.0-2.0*x ;
   tt = t/s2 ;
   mkgegen(z) ;
   for (j=1; j<=maxpol; j++) {  
    q = (double) (j*(j+1))  ; q /= 2.0 ; /* decay */
    xpol = gegen(j, z) ;
    term = ccoeff[j] * xpol * exp (-q*tt) ;
    sum += term ;
    xsum[j-1] = sum ;
    if ((j>=2020) && (term <= 1.0e-10))  {  
     wynn(xsum, j-1, acc, &nacc)   ;
     tt = acc[nacc-1] ;
     if (isnan(tt)) {  
      tt = sum  ;
     }
     if (fabs(tt-sum) < 1.0e-12) break ;
    }
    if (debug) 
    printf("zz %2d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
     j, ccoeff[j], z, xpol, q*tt, term, sum) ;
      if (j==maxpol) { wynn (xsum, maxpol, acc, &nacc)   ;
       tt = acc[nacc-1] ;
       if (isnan(tt)) {  
       tt = sum  ;
      }
     }
   }
   if (first) {  
    first = NO ;
/**
    printmat(xsum, 1, maxpol) ;
    printf("\n\n") ;
    printmat(acc, 1, nacc) ;
*/
   }
// printf("zz %15.9f %15.9f\n", sum, tt) ;
   return tt ;
}

void setccoeff(double p) 
{
   double z, xpol ;
   double top, bot ;
   int j ;

   pco = p ;
   pcmul = p * (1.0-p) ;
   z = 1.0-2*p ;

   for (j=1; j<=maxpol; j++) {  
    top = 4 * (double) (2*j+1) ; 
    bot = (double) (j*(j+1)) ; 
    xpol = gegen(j,z) ; /* efficiency loss here */
    ccoeff[j] = top*xpol/bot ;
//  printf("qqq %9.3f\n", j, ccoeff[j]) ;
   }
}
void mkgegen(double z)  
{
/** A+S Chapter 22  C_n^\alpha(x) */
    int n = maxpol ;
    double alpha = 1.5, x = z ;
    int i ;
    double y0, y1, y2, yn, a, b  ;

    if (z == zco) return ;

    y0 = 1.0 ;
    y1 = 2.0*alpha*x ;
    
    gegenval[0] = y0 ;
    gegenval[1] = y1 ;
    for (i=1; i<n ; i++) {  
     yn = (double) i ;
     a = 2.0*(yn+alpha)*x*y1 ;
     b = (yn+2.0*alpha-1.0)*y0 ;
     y2 = (a-b)/(yn+1.0) ;
     y0 = y1 ; 
     y1 = y2 ;
     gegenval[i+1] = y2 ;
    }
}





double gegen(int n, double x) 
{
/* X_i (z) = T_{i-1} (z) (Kimura) */
 if (x==zco) return gegenval[n-1] ; 
 return cfun(n-1, 1.5, x) ;
}
double cfun (int n, double alpha, double x) 
{
/** A+S Chapter 22  C_n^\alpha(x) */
    int i ;
    double y0, y1, y2, yn, a, b  ;
    y0 = 1.0 ;
    y1 = 2.0*alpha*x ;
    if (n==0) return y0  ;
    if (n==1) return y1  ;
    for (i=1; i<n ; i++) {  
     yn = (double) i ;
     a = 2.0*(yn+alpha)*x*y1 ;
     b = (yn+2.0*alpha-1.0)*y0 ;
     y2 = (a-b)/(yn+1.0) ;
     y0 = y1 ; 
     y1 = y2 ;
    }
    return y2 ;
}

double jac13(int n, double x) 
{
 return jacobi(n, 1, 3, 2*x-1) ;
}


double jacobi (int n, int a , int b, double x) 
{ 
/** A+S Chapter 22  */
    int i ;
    double y0, y1, y2, yn, ya, yb  ;
    double c2, c3, c4 ;
    double cc[3] ;

    y0 = 0.0 ;
    y1 = 1.0 ;
    if (n==0) return y1 ;

    ya = (double) a ;
    yb = (double) b ;

    for (i=0; i<n ; i++) {  
     yn = (double) i ;
     getcc(cc, i, a, b) ;
//   printf("qqq %d %d %d ", i, a, b) ;  printmat(cc, 1, 3) ;
     c2 = cc[0] ;  c3 = cc[1] ; c4 = cc[2]  ;
     y2 = (c2+c3*x)*y1 - c4*y0 ;
     y0 = y1 ; 
     y1 = y2 ;
    }
    return y2 ;
}
double int13(int nn, int ii) 
{
   static double *aco = NULL ;
   double cc[3], dd[3], den, x1[3], x2[3]  ;
   double xm, xn, y ;
   int mp = maxpol+10 ;
   int ns, n, i, j,  k, k0, k1, k2, t  ;
   if (aco==NULL)  {
    ZALLOC(aco, mp*mp, double) ;

    for (ns = 0; ns <= maxpol+2 ; ++ns) {
     for (n=0; n <= ns ; n++) {  
      i = ns -n ;
      if (n>maxpol) continue  ;
      k = kindex(n, i) ;    
      if (n==0) {  
       aco[k] = 1.0 / (double) (i+1) ;
       continue ;
      }
      j = i-4 ;
      if ((i>3) && (j < n)) {   
       k1 = kindex(n, i-1) ;
       aco[k] = aco[k1] ;
       continue ;
      }
      getcc(cc, n-1, 1, 3) ;
      getdd(dd, n-1, 1, 3) ;
      k0 = kindex(n-1,i) ;    
      k1 = kindex(n-1,i+1) ;    
      k2 = kindex(n-2,i) ;    
      x1[0] = aco[k0] ;
      x1[1] = aco[k1] ;
      x1[2] = aco[k2] ;
      dd[2] = -dd[2] ;  
      for (j=0; j<3; j++) {  
       if (x1[j]<0.0) {  
         x1[j] = -x1[j] ;
         dd[j] = -dd[j] ;
       }
      }
      xm = asum(x1,3)/3.0 ;
      vsp(x1, x1, -xm, 3) ;
      y = xm * asum(dd, 3) ;
/** now must add vdot(x1, dd) */
      xn = asum(dd, 3)/ 3.0 ;  
      vsp(dd, dd, -xn, 3) ;
      y += xn * asum(x1, 3) ;
      aco[k]  = y + vdot(x1, dd, 3) ;
/*
      aco[k] = aco[k0]*dd[0]+aco[k1]*dd[1]+aco[k2]*dd[2] ;
*/
     if ( (n<=0) && (i<=3))  {  
/* debug */
       t = (2*n*n+4)*(2*n+2) ; 
       den = (double) t ;
      printf("yy0  %9.3f %9.3f %9.3f\n",cc[0], cc[1], cc[2]) ; 
      printf("zz0 %2d %2d %9.3f %9.3f %9.3f %6d\n",
       n,i,aco[k0], dd[0], dd[0]*den, t ) ; 
      printf("zz1 %2d %2d %9.3f %9.3f %9.3f\n",n,i,aco[k1],dd[1],dd[1]*den) ; 
      printf("zz2 %2d %2d %9.3f %9.3f %9.3f\n",n,i,aco[k2],-dd[2],-dd[2]*den) ; 
     }
     }
    }
    /* end init loop */
   }
    k = kindex(nn, ii) ;    
    return aco[k] ;
}
int kindex(int n, int i) {  
  int m = maxpol+10 ;
  if (n<0) return 0 ; 
  if (i<0) return 0 ; 
  return ((n+1)*m+i) ;
}
void getdd(double *dd, int n, int a, int b) 
{
  /** 
   recursion on 0, 1
  */
  getcc(dd,n, a, b) ;
  dd[0] -= dd[1] ;
  dd[1] *= 2.0 ;
}

void getcc(double *cc, int n, int a, int b) 
{
     double yn,  ya, yb  ;
     double c1, c2, c3, c4 ;

     yn = (double) n ;
     ya = (double) a ;
     yb = (double) b ;

     c1 = 2*(yn+1)*(yn+ya+yb+1)*(2*yn+ya+yb)  ;
     c2 = (2*yn+ya+yb+1)*(ya*ya-yb*yb) ;
     c3 = rfac3(n+n+a+b)  ;
     c4 = 2*(yn+ya)*(yn+yb)*(2*yn+ya+yb+2) ;

     cc[0] = c2 ; cc[1] = c3 ; cc[2] = c4 ;
     vst (cc, cc, 1.0/c1, 3) ;
}
void getcde(double *cc, int n, int a, int b) 
{

/** 
 c, d, e coeffs for 3 -term recurrence
*/
     double yn,  ya, yb  ;
     double c1, c2, c3, c4 ;

     yn = (double) n ;
     ya = (double) a ;
     yb = (double) b ;

     c1 = 2*(yn+1)*(yn+ya+yb+1)*(2*yn+ya+yb)  ;
     c2 = (2*yn+ya+yb+1)*(ya*ya-yb*yb) ;
     c3 = rfac3(n+n+a+b)  ;
     c4 = 2*(yn+ya)*(yn+yb)*(2*yn+ya+yb+2) ;

     cc[0] = c1 ; cc[1] = -c2 ; cc[2] = c4 ;
     vst (cc, cc, 1.0/c3, 3) ;
}
double evalpat(double y, double x, double t, double s2) 
{
  return x*x*(1.0-x)*y*evalpat0(x, y, t, s2) ;
}
double evalpat0(double y, double x, double t, double s2) 
{


   double tt  ;
   double sum = 0.0, term,  q, z, xpol, ypol  ;
   int j, k ;

   tt = t/s2 ;
   for (j=0; j<=maxpol; j++) {  
    k = j+2 ;
    q = (double) (k*(k+1))  ; q /= 2.0 ; /* decay */
    xpol = jac13(j, x) ;
    ypol = jac13(j, y) ;
    term = xpol * ypol * exp (-q*tt) / norm13(j) ;
    if (debug)
     printf("zz %d %12.6f %12.6f %12.6f %12.6f %12.6f\n", 
      j, xpol, ypol, term, sum, norm13(j)) ;
    sum += term ;
   }

   return sum ;
}
double norm11(int n) 
{
   double top, bot ;

   top = (double) (n+1) ;
   bot = (double) (2*n+3)*(n+2) ;

   return top/bot ;
}
double norm13(int n) 
{
/* N_n^13 */
   double top, bot ;

   top = (double) (n+1) ;
   bot = (double) (2*n+5)*(n+4) ;

   return top/bot ;

}

double intx2(double y, int kdeg, double t, double s2) 
{
/** 
 return int (x**kdeg*P_0(x,y,t) d x) 
*/
   double tt  ;
   double sum = 0.0, term,  q, z, xpol, ypol ;
   int j, k ;

   tt = t/s2 ;
   if (tt<=.005) return intx2d(y, kdeg, tt) ;
   for (j=0; j<=maxpol; j++) {  
    k = j+2 ;
    q = (double) (k*(k+1))  ; q /= 2.0 ; /* decay */
    xpol = int13(j, kdeg) ;
    ypol = jac13(j, y) ;
    term = xpol * ypol * exp (-q*tt) / norm13(j) ;
    sum += term ;
   }
   return sum ;
}
double intx2d(double y, int kdeg, double tt) 
{
/** 
 int x**kdeg*C(x,y; t) / y*x*x*(1-x)) ;
 differential approximation for small tt
*/
 double top, bot, valbase, t1, t2, der, z1 ;

 if ((y==0.0) || (y==1.0)) fatalx("(intx2d) bad y\n") ;
  top = pow(y, (double) kdeg);
  bot = y*y*y*(1.0-y) ;
  valbase = top/bot ;
  z1 = 0.5 * (double) (kdeg-1)*(kdeg-2)  ;
  t1 = z1*pow(y, (double) (kdeg-3) )/y ;   
  t2 = valbase/y ;

  der = t1 -t2 ; 
  return valbase + tt*der ;
}


double intx(double y, int kdeg, double t, double s2) 
{
/** 
 return int (x**kdeg*P_0(y,t,x) d x) 
*/
   double tt  ;
   double sum = 0.0, term,  q, z, xpol  ;
   int j, k ;

   if (y != jpco) {  
    setjcoeff(y) ;
   }
   tt = t/s2 ;
   for (j=0; j<=maxpol; j++) {  
    k = j+2 ;
    q = (double) (k*(k+1))  ; q /= 2.0 ; /* decay */
    xpol = int13(j, kdeg) ;
    term = jcoeff[j] * xpol * exp (-q*tt) ;
    if (term>100000)  {  
     printf("badz %d %d %9.3f %g %g %g\n",kdeg, j, tt, xpol, jcoeff[j], exp (-q*tt) ) ;
    }
    sum += term ;
   }
   return sum ;
}

void setjcoeff(double p) 
{
/** 
 this is tricky  
 we define P(x,p,t) to be backward transition prob.  
 evaluated by evalpat  
 and P_0(x,p,t) = (x,p,t)/(p*p*(1-p))

 we compute jcoeff as appropriate expansion of P_0 and an 
 extra factor of p gets in since the weight funcion has an extra 
 factor of p 
*/
 
   double  xpol, x ;
   double top, bot ;
   int j, k ;

   jpco = p ;
   jpcmul = p *  p * (1.0-p) ;

   x = p ;
   for (j=0; j<=maxpol; j++) {  
    bot =  (double) (j+1) ; 
    k = (j+j+5)*(j+4) ;  
    top = p * (double) (k) ; 
    xpol = jac13(j,x) ; /* efficiency loss here */
    jcoeff[j] = top*xpol/bot ;
    if (j<=3) printf("zzjco %d %9.3f %12.6f\n", j, p, jcoeff[j]) ;
   }
}
double rfac3 (int m) 
{
 int mm ;
 mm = m*(m+1)*(m+2) ;
 return (double) mm ;
}
double probfix(double y, double t) 
/* 
 return probability allele is extinct by time t
 probfix(1-y, t) is prob allele fixes
*/
{
  return 0.5*y*(1.0-y)*kimint0(y,t,NO);
}

double kimint0(double y,  double t, int isupper)               
{

/** 
 \int_t^infinity K_0 (0, y, s) d s
 for isupper = YES 

 \int_0^t   K_0 (0, y, s)  for isupper = NO

 
*/

   double tt, s2= 1.0, x=0.0, p=0.0   ;
   double sum = 0.0, term, alpha, q, z, xpol  ;
   int j ;

   if (isupper == NO) {  
    sum = (2.0/y) - kimint0(y, t, YES) ;
    return sum ;
   }
   if (t < .000001) return 2.0/y ; 
/* should use a different expansion really */
   if (p != pco) {  
    setccoeff(p) ;
   }
   z = 1.0-2.0*y ;
   tt = t/s2 ;
   mkgegen(z) ;
   for (j=1; j<=maxpol; j++) {  
    q = (double) (j*(j+1))  ; q /= 2.0 ; /* decay */
    xpol = gegen(j, z) ;
    term = ccoeff[j] * xpol * exp (-q*tt) ;
    sum += term/q ;
   }

   return sum ;
}

void setmaxpol(int maxp) { 
 if ((maxp <= 0) || (maxp>MAXPOL)) maxpol = MAXPOL ;
 else maxpol = maxp ;

}

double getccoeff(int j) 
{
  return ccoeff[j] ;
}

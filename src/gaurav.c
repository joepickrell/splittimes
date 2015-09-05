#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  

#include <nicklib.h>  
#include "kimsubs.h" 

int verbose = NO ;
char *iname = "gdata" ;
void readcommands(int argc, char **argv) ;

void doit(double fst, double p1, double p2) ;


int main  (int argc , char **argv) 
{
   double x, t, y, ktrans, tau ;
   int i , j ;
   int n, nx, kx, ndata ; 
   double **xx ; 
   double y0, yn, z ;
   double p0, p1, u1, u2, urat, pp0 ;
   double *udata, *xval, *ww ;

  readcommands(argc, argv) ;
  n = numlines(iname) ; 
  xx = initarray_2Ddouble(3, n, 0.0) ;

  n =  getxx(xx, n, 3, iname) ;

  for (i=0; i<n ; ++i) { 
   doit(xx[0][i], xx[1][i], xx[2][i]) ;
  }
  printf("## end of run\n") ;
  return 0 ;
}

void readcommands(int argc, char **argv) 

{
  int i ;

  while ((i = getopt (argc, argv, "i:V")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

}

void doit(double fst, double p1, double p2) 
{
  double tau ;
  double var, sig, kimd, gaussd ;
  static int ncall = 0 ;

  ++ncall ;

  if (ncall==1) { 
   printf("%12s", "Fst") ;
   printf("  %12s", "p1") ;
   printf("  %12s", "p2") ;
   printf(  "%12s", "Gaussian") ;
   printf(  "%12s", "Kimura") ;
   printnl() ;
  }

  printf("%12.6f", fst) ;
  printf("  %12.6f", p1) ;
  printf("  %12.6f", p2) ;
  tau = -log(1-2*fst) ;  // fst < 1/2
  kimd = evalkim(p2, p1, tau, 1.0) ;
  var = 2.0*fst*p1*(1.0-p1) ;
  sig = sqrt(var) ;
  gaussd = nordis( (p2-p1)/sig ) /sig ;
  printf(  "%12.6f", gaussd) ;
  printf(  "%12.6f", kimd) ;
  printnl() ;
 
}


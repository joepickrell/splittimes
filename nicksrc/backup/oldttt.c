#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <nicklib.h>

void twt(double x)  ;
double twtail(double twstat)  ;

main()
{
  twt(0.1) ;
  twt(-1.2) ;
  twt(4.8) ;
  twt(-12.0) ;
  twt(12.0) ;
}
void twt(double x) 
{
  double y ;
  y = twtail(x) ;
  printf ("%12.6f  %12.6f\n", x, y) ;

}
double twtail(double twstat) 
{
#define TWTAB  "/home/np29/tables/twtable.txt"  
#define TABSIZE 1000 
#define MAXSTR  128 
 static double inc = -1 ;
 static double xval[TABSIZE], xpdf[TABSIZE], xtail[TABSIZE] ;
 char sa[40], sb[40] ;
 double tail, y1, y2 ; 
 FILE *fff ;
 char line[MAXSTR] ;
 static int num=0 ;
 int k ;


  if (num == 0)  {
  fff = fopen(TWTAB, "r") ;
  if (fff==NULL) fatalx("can't open file %s \n",TWTAB) ;
  
  while (fgets(line, MAXSTR, fff) != NULL)  { 
    sscanf(line, "%s %s", sa, sb) ;
    xval[num] =  atof(sa) ;
    xpdf[num] = atof(sb) ;
    ++num ;
   }
   inc = xval[1] - xval[0]   ;
   fclose(fff) ;
   printf("twtab loaded: %d  %12.6f\n", num, inc ) ;
   tail = 0.0 ;
   for (k=num-1; k>0; --k) {  
     xtail[k] = tail ;
     tail += inc*0.5*(xpdf[k]+xpdf[k-1]) ;
   }
   xtail[0] = 1.0 ;
   printf ("##check: %12.6f\n", xtail[1]) ;
  }
  k = firstgt(twstat, xval, num) ;
  if (k==0) return 1.0 ;
  if (k>=num) return 0.0 ;
  y1 = (xval[k]-twstat)/inc ;
  y2 = 1.0-y1 ;
  return y1*xtail[k-1] + y2*xtail[k] ;
}


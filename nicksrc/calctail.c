#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  

#include <nicklib.h> 

#define MAXFIELD 10 
#define MAXS  512

int main() 
{
  double p, val ;
  int n, t, c ;
  char str[MAXS] ;
  char line[MAXS] ;
  char *spt[MAXFIELD] ;
  int nsplit ;

  while (fgets(line,MAXS,stdin) != NULL)   {
   nsplit = splitup(line, spt, MAXFIELD) ;
   if ((nsplit <3) || (nsplit>4)) fatalx("bad line %s\n",line) ;
   c = '+' ;
   if (nsplit == 4) c = spt[3][0] ;
   n = atoi(spt[0]) ;
   t = atoi(spt[1]) ;
   p = atof(spt[2]) ;
   freeup(spt, nsplit) ;
   if ((p<=0.0) || (p>=1.0)) fatalx("bad line %s\n",line) ;
   val = binomtail(n,t,p,c) ;
   printf ("%15.6e\n",val) ;
  }
}

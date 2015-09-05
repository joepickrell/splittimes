#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <stdarg.h>
#include <sys/times.h>
#include <time.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>


#define MAXSTR 10000
#define MAXFF 50

#include "strsubs.h" 
#include "vsubs.h" 

extern int errno ;

int splitup(char *strin, char**spt, int maxpt) 
/**
 improved code
 oldsplitup for old code.  
 retained in case there are compatibility problems 
*/
{
 char *s1, *s2, *sx ;
 char *str ;
 int i, len, num ;

 len = strlen(strin) ;
 if (len==0) return 0 ;  
 ZALLOC(str, 2*len, char) ; 
 for (num=0,  i=0, sx=strin; i<maxpt; i++)  {
  s1 = fnwhite(sx) ;
  if (s1==NULL) { 
   break ;
  }
  s2 = fwhite(s1) ;
  if (s2==NULL) {
   s2 = s1+strlen(s1) ;
  }
  s2-- ;  /* now points at last character of next word */
  len = s2-s1+1 ;
  strncpy(str,s1,len) ;
  str[len] = '\0' ;
  spt[num] = strdup(str) ;
  ++num ;
  sx = s2+1 ;
 }
 free(str) ;
 return num ;
}

void freeup (char *strpt[],int numpt)
/** free up array of strings */
{
    int i ;
    for (i=0; i<numpt; i++)  {
      if (strpt[i] != NULL) freestring(&strpt[i]) ;
    }
}

int first_word(char *string, char *xword, char *xrest)

/*  first_word(string, *word, *rest)

        Break the string into the first word and the rest.  Both word and
rest begin with non-white space, unless rest is null.
        Return:
                0 means string is all white
                1 means word is non-white, but rest is white
                2 means word and rest are non-white
 
        If string and rest coincide, string will be overwritten
 
*/
{
   char *spt, x ;
   char *ss = NULL, *sx ;
   int l1,l2 ;

   ss = strdup(string) ;
   if (ss==NULL)  {
    printf("strdup fails\n") ;
    printf("%s\n",string) ;
    fatalx("first_word... strdup fails\n") ;
   }
   fflush(stdout) ;
   spt = ss ;
   xword[0]=xrest[0] = '\0' ;
   if ((spt = fnwhite(ss)) == NULL) {
     free(ss) ;
     return 0 ;
   }
    sx = fwhite(spt) ;
    if (sx==NULL) {
     strcpy(xword,spt) ;
     free (ss) ;
     return 1 ;
    }
    l1 = sx-spt ;
    l2 = strlen(sx) - 1 ;
    *sx = '\0' ;
    strcpy(xword,spt) ;
    if (l2 <= 0) {
     free(ss) ;
     return 1 ;
    }

    sx = fnwhite(sx+1) ;
    if (sx==NULL) { 
     free (ss) ;
     return 1 ;
    }
    strcpy(xrest,sx) ;  
    free (ss) ;
    return 2 ;
}

char * fnwhite (char * ss) 
/* return first non white space */
{
  char *x ;
  if (ss==NULL) fatalx("fnwhite: logic bug\n") ;
  for (x= ss; *x != '\0'; ++x)  {
   if (!isspace(*x)) return x ;
  }
  return NULL ;
}

char * fwhite (char * ss) 
/* return first white space */
{
  char *x ;
  int n ;
  for (x= ss; *x != '\0'; ++x)  {
   if (isspace(*x)) return x ;
  }
  return NULL ;
}
 
static char Estr[MAXSTR];

void fatalx( char *fmt, ...) 
{	va_list args;

	va_start( args, fmt);
	vsprintf( Estr, fmt, args);
	va_end( args);
        fflush(stdout) ;

        fprintf(stderr,"fatalx:\n%s",Estr) ;
        fflush(stderr) ;
        abort() ;
}
int NPisnumber (char c) 
/**
 returns 1 if - + or digit 
*/
{
  if (isdigit(c)) return 1 ;
  if (c=='+') return 1 ;
  if (c=='-') return 1 ;

  return 0 ;
}
int isnumword (char *str) 
{

  int i, len, numpt ;
  char c ;
  len = strlen(str) ;
 
  numpt = 0 ;
  for (i=0; i<len; i++) {
   c = str[i] ; 

   if ((c == '.') && (numpt==0) ) {
    ++numpt ;
    continue ;
   }

   if (!NPisnumber(c)) return NO ;
   if (!isdigit(c) && (i>0) ) return NO ;

  }
  return YES ;

}

long seednum() 
{
   long a, b, c, d ;
   struct tms tbuff ;

   a = (long) getpid() ;
   b = (long) getuid() ;
   d = times(&tbuff) ;

   c = d ^ ((a + b) << 15) ;


   return c ;

}
int split1 (char *strin,  char *strpt[], char splitc)
/*
take a string and break it into 2 substrings separated by splitc ;
\
numpt is number of words returned  (1 or 2) 
*/
{
  char rest[MAXSTR],str[MAXSTR],ww[MAXSTR]  ;
  int  len, i, l ;

  strncpy(str,strin,MAXSTR) ; 
  len = strlen(strin) ;
  for (i=0; i<len; i++) {
   if (str[i] == splitc) {
    l = i ;
    strncpy(ww, str, l) ;  
    ww[l] = '\0' ;
    strpt[0] = strdup(ww) ;
    l = len-(i+1) ;
    if (l<=0) return 1 ;
    strncpy(rest,str+i+1,l) ;
    rest[l] = '\0' ;
    strpt[1] = strdup(rest) ;
    return 2 ;
   }
  }
  strpt[0] = strdup(strin) ;
  return 1 ;
}

int oldsplitup (char *strin,  char *strpt[],int maxpt)
/*
take a string and break it into words;      \
numpt is number of words returned max maxpt
if numpt==maxpt then last word may still need
splitting
*/
{
  char rest[MAXSTR],str[MAXSTR],ww[MAXSTR]  ;
  int  numpt ;

  strncpy(str,strin,MAXSTR) ; numpt=0 ;
  while (first_word(str,ww,rest)>0) {
  strpt[numpt] = strdup(ww) ;
  strncpy(str,rest,MAXSTR) ;
  numpt++ ;

  if (numpt==maxpt)  break ;
  }
  return numpt ;

}

void printbl(int n) 
{
  int i ;
  for (i=0; i<n; i++) {
   printf(" ") ;
  }
}

void printnl() 
{
   printf("\n") ;
}

void striptrail(char *sss, char c) 
/** 
 strip out trailing characters 
 c will usually be ' '
*/
{
   int len, i ;
   len = strlen(sss) ;
   for (i=len-1; i>=0; --i) {  
    if (sss[i] != c) return ;
    sss[i] = '\0' ;
   }
}

void catx(char *sxout, char **spt, int n) 
{
     int i  ;
     sxout[0] = CNULL ; 

     for (i=0; i<n; i++) {
      strcat(sxout,spt[i]) ;
     }


}

void catxx(char *sxout, char **spt, int n) 
/** 
 like catx but with space between items 
*/
{
     int i  ;
     sxout[0] = '\0' ; 

     for (i=0; i<n; i++) {
      strcat(sxout,spt[i]) ;
      if (i<(n-1)) strcat(sxout, " ") ;
     }
}

void makedfn(char *dirname, char *fname, char *outname, int maxstr) 
/** makes full path name.    
  If fname starts with '/' or dirname = NULL we 
  so nothing. 
  outname MUST be allocated of length at least maxstr 
*/
{
    char *ss ;
    int len ;
 
    if ((dirname==NULL) || (fname[0]=='/')) {  
/* if fname starts with / we assume absolute pathname */
     len = strlen(fname) ;  
     if (len>=maxstr) fatalx("(makedfn) maxstr too short\n") ;
     strcpy(outname, fname) ;
     return ;
    }
    len = strlen(dirname) + strlen(fname)+1 ;  
    if (len>=maxstr) fatalx("(makedfn) maxstr too short\n") ;

    ss = outname ;
    strcpy(ss,dirname) ;
    ss = ss+strlen(dirname) ;  
    ss[0] = '/' ;
    ++ss ;
    strcpy(ss,fname) ;
}

int substringx (char **ap, char *inx, char *outx, int niter)
/** 
 *ap is original string 
 all occurrences of inx are substituted with outx 
 can loop so be careful !!  

 NB.  ap must be on heap.  Fixed allocation not supported 
*/
{
   char *a, *pt ;
   char *str ;
   int len, off, x ;

   if (niter>50) fatalx("bad string replacement\n %s\n", *ap) ;

   a = *ap ;
   len = strlen(a) + strlen(inx) + strlen(outx) + 1 ;
   pt = strstr(a, inx) ;
   if (pt == NULL) { 
    return 0 ;
   }
   ZALLOC(str, len, char) ;
   off = pt - a ;
   strncpy(str,a,off) ;            
   strcpy(str+off, outx) ;
   x = strlen(outx) ; 
   pt += strlen(inx) ;
   strcpy(str+off+x, pt) ;
   
   free(a) ;
   *ap = strdup(str) ;
   free(str) ;
   return (1 + substringx(ap, inx, outx, niter+1)) ;
}

int substring (char **ap, char *inx, char *outx)
/** 
 *ap is original string 
 all occurrences of inx are substituted with outx 
 can loop so be careful !!  
 

 NB.  ap must be on heap.  Fixed allocation not supported 
*/
{
   return (substringx(ap, inx, outx, 0)) ;
}



int numcols(char *name)
// number of cols 
{
  FILE *fff ;
  char line[MAXSTR] ;
  char *spt[MAXSTR] ;
  char *sx ;
  int nsplit, num=0 ;

  if (name == NULL) fatalx("(numlines)  no name")  ;
  openit(name, &fff, "r") ;
  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, MAXFF) ; 
   if (nsplit==0) continue ;
   sx = spt[0] ;
   if (sx[0] == '#') {
    freeup(spt, nsplit) ;
    continue ;
   }
   freeup(spt, nsplit) ;
   fclose(fff) ;
   return nsplit ;
  }
}

int numlines(char *name)
// number of lines   no comments or blanks
{
  FILE *fff ;
  char line[MAXSTR] ;
  char *spt[MAXSTR] ;
  char *sx ;
  int nsplit, num=0 ;

  num = 0; 
  if (name == NULL) fatalx("(numlines)  no name")  ;
  openit(name, &fff, "r") ;
  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, MAXFF) ; 
   if (nsplit==0) continue ;
   sx = spt[0] ;
   if (sx[0] == '#') {
    freeup(spt, nsplit) ;
    continue ;
   }
   ++num ;
   freeup(spt, nsplit) ;
  }
  fclose(fff) ;
  return num ;
}

void openit(char *name, FILE **fff, char *type)  
{
  char *ss ;  
  if (name==NULL) fatalx("(openit) null name\n") ;
  *fff = fopen(name,type) ;
  if (*fff==NULL) {
   ss = strerror(errno) ;
   printf("bad open %s\n, name") ; 
// system("lsof | fgrep np29") ;
   fatalx("can't open file %s of type %s\n error info: %s\n",name,type,ss) ;
  }
}

int 
getxx(double **xx, int maxrow, int numcol, char *fname)
{

  char line[MAXSTR] ;
  char *spt[MAXFF] ;
  char *sx ;
  int nsplit, i, j, num=0, maxff ;
  FILE *fff ;
  int nbad = 0 ; 

  if (fname == NULL) fff = stdin ; 
  else {
   openit(fname, &fff, "r") ;
  }
  maxff = MAX(MAXFF, numcol) ; 

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, maxff) ; 
   sx = spt[0] ;
   if (sx[0] == '#') {
    freeup(spt, nsplit) ;
    continue ;
   }
   if (nsplit<numcol) { 
     ++nbad ;
     if (nbad<10) printf("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol, line) ;
     continue ;
   }
   if (num>=maxrow) fatalx("too much data\n") ;
   for (i=0; i<numcol; i++)  {
    xx[i][num]  = atof(spt[i]) ;
   }
   freeup(spt, nsplit) ;
   ++num ;
  }
  if (fname != NULL) fclose(fff) ;
  return num ;
}

double clocktime()
{
  double xtime ;
  double y ;

  xtime = (double) clock() ;
  y =  xtime / (double) CLOCKS_PER_SEC ;
  return y ;
}
int
indxstring(char **namelist, int len, char *strid)  
// look for string in list.  Was called indxindex
{
     int k ; 
     for (k=0; k< len; k++) {  
      if (strcmp(namelist[k], strid) == 0) return k ;
     }
     return -1 ;
}

int 
getxxnames(char ***pnames, double **xx, int maxrow, int numcol, char *fname)
{

  char line[MAXSTR] ;
  char *spt[MAXFF] ;
  char *sx ;
  int nsplit, i, j, num=0, maxff, numcolp ;
  FILE *fff ;
  int nbad = 0 ; 
  char **names ;

  names = *pnames ;
  if (fname == NULL) fff = stdin ; 
  else {
   openit(fname, &fff, "r") ;
  }
  numcolp = numcol + 1 ;
  maxff = MAX(MAXFF, numcolp) ; 

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, maxff) ; 
   if (nsplit == 0) { 
    freeup(spt, nsplit) ;
    continue ;
   }
   sx = spt[0] ;
   if (sx[0] == '#') {
    freeup(spt, nsplit) ;
    continue ;
   }
   names[num] = strdup(sx) ;
   if (nsplit<numcolp) { 
     ++nbad ;
     if (nbad<10) printf("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol, line) ;
     continue ;
   }
   if (num>=maxrow) fatalx("too much data\n") ;
   for (i=0; i<numcol; i++)  {
    xx[i][num]  = atof(spt[i+1]) ;
   }
   freeup(spt, nsplit) ;
   ++num ;
  }
  if (fname != NULL) fclose(fff) ;
  return num ;
}

int 
getxxnamesf(char ***pnames, double **xx, int maxrow, int numcol, FILE *fff)
/** 
like getxxnames but file already open 
*/
{

#define MAXFF  50

  char line[MAXSTR] ;
  char *spt[MAXFF] ;
  char *sx ;
  int nsplit, i, j, num=0, maxff, numcolp ;
  int nbad = 0 ; 
  char **names ;

  names = *pnames ;

  numcolp = numcol + 1 ;
  maxff = MAX(MAXFF, numcolp) ; 

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, maxff) ; 
   if (nsplit == 0) { 
    freeup(spt, nsplit) ;
    continue ;
   }
   sx = spt[0] ;
   if (sx[0] == '#') {
    freeup(spt, nsplit) ;
    continue ;
   }
   names[num] = strdup(sx) ;
   if (nsplit<numcolp) { 
     ++nbad ;
     if (nbad<10) printf("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol, line) ;
     continue ;
   }
   if (num>=maxrow) fatalx("too much data\n") ;
   for (i=0; i<numcol; i++)  {
    xx[i][num]  = atof(spt[i+1]) ;
   }
   freeup(spt, nsplit) ;
   ++num ;
  }
  return num ;
}


int 
getss(char **ss, char *fname)
/** 
 get list of names 
*/
{

  char line[MAXSTR] ;
  char qqq[MAXSTR] ;
  char *spt[MAXFF] ;
  char *sx ;
  int nsplit, i, j, num=0, maxff ;
  FILE *fff ;


  if (fname == NULL) fff = stdin ; 
  else {
   openit(fname, &fff, "r") ;
  }
  maxff = MAXFF  ; 

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, maxff) ; 
   sx = spt[0] ;
   if (sx[0] == '#') {
    freeup(spt, nsplit) ;
    continue ;
   }
   if (nsplit<1) { 
     continue ;
   }
   ss[num] = strdup(spt[0]) ;
   freeup(spt, nsplit) ;
   ++num ;
  }
  if (fname != NULL) fclose(fff) ;
  return num ;
}

char *int_string(int a, int len, int base)
{
  static char ss[100] ;
  int t = a, k, i ;
  char *binary = "01" ;

  ss[len] = CNULL ;
  for (i=0; i<len; i++) {
   k = t % base  ;
   ss[len-i-1] = '0' + k ;   
   t = t/base ;
  }
  return ss ;
// fragile
}
char *binary_string(int a, int len)
{
  static char ss[100] ;
  int t = a, k, i ;
  char *binary = "01" ;

  ss[len] = CNULL ;
  for (i=0; i<len; i++) {
   k = t % 2 ;
   ss[len-i-1] = binary[k] ;
   t = t/2 ;
  }
  return ss ;
// fragile
}
char num2base (int num) 
{

  char *bases ="ACGT", c ;
  c = '?' ;
  if (num<0) return c ;
  if (num>3) return c ;
  return bases[num] ;

}

int base2num(char c)

{
    char cc ; 

    cc = toupper(c) ;

    switch (cc)  {
     case 'A':  return 0;
      break ;
     case 'C':  return 1;
      break ;
     case 'G':  return 2;
      break ;
     case 'T':  return 3;
      break ;
     default:  return -1 ;
    }
}

int string_binary(char *sx)
{
   int *aa, len, i, t ;
   char c ;

   len = strlen(sx) ;
   ZALLOC(aa, len, int) ;

   for (i=0; i<len; i++) {

    c = sx[i] ;
    if (c == '0') continue ;
    if (c != '1') fatalx("bad string: %s\n", sx) ;
    aa[i] = 1;
   }
   t = kodeitb(aa, len, 2) ;
   free(aa) ;
   return t ;

}
void freestring(char **ss)  
/* note extra indirection */
{
  if (*ss == NULL) return ; 
  free(*ss) ; 
  *ss = NULL ;
}

void copystrings(char **sa, char **sb, int n) 
{
  int i ; 
  for (i=0; i<n; ++i) {  
   sb[i] = strdup(sa[i]) ;
  }
}
void printstrings(char **ss, int n) 
{
  int k ;

  for (k=0; k<n; ++k) {
   printf("%s\n", ss[k]) ;
  }
}

int ridfile(char *fname)
{
  int t ;

  chmod(fname, 0777) ;
  t = unlink(fname) ;
  return t ;
}

char compbase(char x) 
// upper case !!
{
 if (x=='A') return 'T' ;
 if (x=='C') return 'G' ;
 if (x=='G') return 'C' ;
 if (x=='T') return 'A' ;

 return x ;

}



int 
getjjnames(char ***pnames, int **jj, int maxrow, int numcol, char *fname)
{

  char line[MAXSTR] ;
  char *spt[MAXFF] ;
  char *sx ;
  int nsplit, i, j, num=0, maxff, numcolp ;
  FILE *fff ;
  int nbad = 0 ; 
  char **names ;

  names = *pnames ;
  if (fname == NULL) fff = stdin ; 
  else {
   openit(fname, &fff, "r") ;
  }
  numcolp = numcol + 1 ;
  maxff = MAX(MAXFF, numcolp) ; 

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, maxff) ; 
   if (nsplit == 0) { 
    freeup(spt, nsplit) ;
    continue ;
   }
   sx = spt[0] ;
   if (sx[0] == '#') {
    freeup(spt, nsplit) ;
    continue ;
   }
   names[num] = strdup(sx) ;
   if (nsplit<numcolp) { 
     ++nbad ;
     if (nbad<10) printf("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol, line) ;
     continue ;
   }
   if (num>=maxrow) fatalx("too much data\n") ;
   for (i=0; i<numcol; i++)  {
    jj[i][num]  = atoi(spt[i+1]) ;
   }
   freeup(spt, nsplit) ;
   ++num ;
  }
  if (fname != NULL) fclose(fff) ;
  return num ;
}

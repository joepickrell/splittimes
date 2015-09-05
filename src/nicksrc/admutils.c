#include "admutils.h" 

void loadstats(FILE *statsfile, Indiv *indiv_array, int *numloaded)
{
  int i,c, numcounted, numl=0;
  
  numcounted = countlines(statsfile)-1;
  while ( (c = getc(statsfile)) != '\n' ) {}
  
  for (i=0; i<numcounted; i++) {
    /* printf("%d\t%d\n",i,countcolumns(statsfile)); */
    indiv_array[i].ignore = NO ;
    if (countcolumns(statsfile) == 6) {
      if ( fscanf(statsfile, "%d %s %1s %lf %lf %*lf", 
		  &(indiv_array[i].affstatus), indiv_array[i].ID,
		  &(indiv_array[i].gender),
		  &(indiv_array[i].theta_mode), &(indiv_array[i].lambda_mode)) == 5 ) {
	numl++;
      }
    }
    else if (countcolumns(statsfile) == 11) {
      if ( fscanf(statsfile, "%d %s %1s %lf %lf %lf %lf %*lf %lf %lf %*lf", 
		  &(indiv_array[i].affstatus), indiv_array[i].ID,
		  &(indiv_array[i].gender),
		  &(indiv_array[i].theta_mean), &(indiv_array[i].theta_stderr),
		  &(indiv_array[i].lambda_mean), &(indiv_array[i].lambda_stderr),
		  &(indiv_array[i].theta_mode), &(indiv_array[i].lambda_mode)) == 9 ) {
	numl++;
      }
    }
    else {
      printf("loadstats error with individual %d\n",i); 
      continue;
    }
    while ( (c = getc(statsfile)) != '\n' ) {}
  }
  if (numcounted != numl) {
    fatalx("bad statsfile numcounted: %d numloaded: %d\n",numcounted,numl);
  }
  *numloaded = numl ;
}

void loadXstats(FILE *Xstatsfile, Indiv *indiv_array, int numindivs, int *numloaded)
     /* numindivs is size of indiv_array;
      *numloaded is number of individuals loaded from Xstatsfile */
{
  int i,j,c, numcounted, numcols, affstatus=0;
  int numl = 0 ;
  char ID[30];

  numcounted = countlines(Xstatsfile)-1;
  while ( (c = getc(Xstatsfile)) != '\n' ) {}
  
  numcols = countcolumns(Xstatsfile);

  for (i=0; i<numcounted; i++) {
    fscanf(Xstatsfile, "%d %s", &affstatus, ID);
    for (j=0; j<numindivs; j++) {
      if (strcmp(indiv_array[j].ID,ID)==0) break;
    }
    if (j>=numindivs) fatalx("(loadX) indiv %s not found\n", ID) ;
    indiv_array[j].affstatus = affstatus;
    if (numcols == 6) {
      if ( fscanf(Xstatsfile, "%1s %lf %lf %*lf", 
		  &(indiv_array[j].gender),
		  &(indiv_array[j].Xtheta_mode), &(indiv_array[j].Xlambda_mode)) == 3 ) {
	numl++;
      }
    }
    else if (numcols == 11) {
      if ( fscanf(Xstatsfile, "%1s %lf %lf %lf %lf %*lf %lf %lf %*lf", 
		  &(indiv_array[j].gender),
		  &(indiv_array[j].Xtheta_mean), &(indiv_array[j].Xtheta_stderr),
		  &(indiv_array[j].Xlambda_mean), &(indiv_array[j].Xlambda_stderr),
		  &(indiv_array[j].Xtheta_mode), &(indiv_array[j].Xlambda_mode)) == 7 ) {
	numl++;
      }
    }
    while ( (c = getc(Xstatsfile)) != '\n' ) {}
  }
  if (numcounted != numl) {
    fatalx
     ("badXstatsfile numcounted: %d\t numloaded: %d\n",numcounted,numl);
  }
  *numloaded = numl ;
}

int countlin(char *fname) {
 FILE *fp ;
 int t ;

 openit(fname, &fp, "r") ;
 t = countlines(fp) ;
 fclose(fp) ;
 return t ;
}

int countcol(char *fname) {
 FILE *fp ;
 int t ;

 openit(fname, &fp, "r") ;
 t = countcolumns(fp) ;
 fclose(fp) ;
 return t ;
}


int countlines(FILE *fp) 
{     /* count number of newlines followed by a graphic character */
  int i=1,c,linetest;

  rewind(fp); 
  while ( (c = getc(fp)) != EOF ) {
    if (c == '\n') {
      linetest = 0;
      while (c = getc(fp)) {
	if (c == '\n' || c == EOF) break;
	if (isgraph(c)) { linetest=1; break; }
      }
      if (linetest) i++;
      else break;
    }
  }
  rewind(fp);
  return i;
}

int countcolumns(FILE *fp)
{  /* count number of text columns separated by whitespace */
  int i=0,c;
  fpos_t ptr;

  if (fgetpos(fp,&ptr)) {
    printf("error counting columns\n");
    return 0;
  }

  while ( (c = getc(fp)) != '\n' ) {
    if (isgraph(c)) {
      i++;
      while (isgraph(c = getc(fp))) {}
      ungetc(c,fp);
    }
  }
  fsetpos(fp,&ptr);
  return i;
}

void sett1(double *tt, double theta, int numstates) 
{
    if (numstates==2)  { 
     tt[0] = 1.0-theta ;
     tt[1] = theta ;
     tt[2] = 0.0 ;
    }
    if (numstates==3)  { 
     tt[0] = (1.0-theta)*(1.0-theta) ;
     tt[1] = 2.0*theta*(1.0-theta) ;
     tt[2] = theta*theta ;
    }
}

void sett1r(double *tt, double theta, int numstates, double risk)
{
  double y ;
  sett1(tt, theta, numstates) ;
  tt[1] *= risk ;
  tt[2] *= risk*risk ;  
  y = asum(tt, numstates) ;
  vst(tt, tt, 1.0/y, numstates) ;
}

void gettln(SNP *cupt, Indiv *indx, 
  double *ptheta, double *plambda, int *pnumstates, int *pignore) 
/* set theta, lambda numstates */
{
   double theta, lambda ; 
   int numstates, chrom, ignore ;

   theta  = indx->theta_mode ;
   lambda = indx->lambda_mode ;
   ignore = indx->ignore ;
   numstates = 3 ;

   chrom = cupt -> chrom ;

   if (chrom == 23) { 
    theta  = indx->Xtheta_mode ;
    lambda = indx->Xlambda_mode ;
    ignore = indx->Xignore ;
    if (indx -> gender == 'M') numstates = 2;
   }
   *ptheta = theta ;
   *plambda = lambda ;
   *pnumstates = numstates ;
   if (pignore != NULL) 
    *pignore = ignore ;
}


void puttln(SNP *cupt, Indiv *indx, 
  double ptheta, double plambda) 
/* put theta, lambda */
{

   int chrom ;


   chrom = cupt -> chrom ;

   if (chrom == 23) { 
    indx->Xtheta_mode = ptheta;
    indx->Xlambda_mode  = plambda ;
    return ;
   }
   indx->theta_mode = ptheta;
   indx->lambda_mode  = plambda ;
   return ;
}


/******** UTILITY FUNCTIONS **********/
double **initarray_2Ddouble(int numrows, int numcolumns, double initval)
{
  int i,j;
  double **array;

   		
  ZALLOC(array, numrows, double *) ;
  for (i=0; i<numrows; i++) {
    ZALLOC(array[i],numcolumns,double);
    if (initval != 0.0) 
     vclear(array[i], initval, numcolumns) ;
  }
  return array;
}
void free2D (double ***xx, int numrows) 
{
   double **array ;
   int i ;
   array = *xx ;

  for (i=0; i<numrows; i++) {
    free(array[i]) ;
  }
  free(array) ;
}


void fataly(const char *name)
{
  printf("%s",name);
  exit(1);
}

int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

void openit(char *name, FILE **fff, char *type)  
{
  if (name==NULL) fatalx("(openit): null name\n") ;
  *fff = fopen(name,type) ;
  if (*fff==NULL) fatalx("can't open file %s of type %s\n",name,type) ;
}

void pcheck (char *name, char x) 
{
 
 if (name == NULL) 
  if (x != NULL) 
   fatalx ("parameter %c compulsory\n",x) ;
  else fatalx("missing argument\n") ;
}

void printm(double **M, int numstates)  
{
     int i,j ;
     printf("M:\n") ;
     for (i=0; i<numstates; i++)  {  
      for (j=0; j<numstates; j++)  {  
       printf("%9.3f ", M[j][i]) ;
      }
      printf("\n") ;
     }
}


int numvalidgtypesx(int indiv, SNP **snpm,  int fc, int lc)
/* count valid gtypes for individual indiv */
{
  int nvalid,k, n ;
  SNP *cupt ;
  nvalid = 0 ;
  for (k=fc; k <= lc; k++) {
   cupt = snpm[k] ;
   if (cupt -> isfake)  continue ;
   if (cupt -> ignore) continue ;
   n = cupt -> ngtypes ;
   if (n==0) continue ;
   if (indiv >= n) fatalx("(numvalidgtypesx) bad index %d %d\n",
     k, indiv) ;
   if (cupt->gtypes[indiv] >= 0) ++nvalid ;
  }
  return nvalid ;
} 

int numvalidgtypes(SNP *cupt) 
{
  int nvalid, n, i, k ;
  if (cupt -> isfake) return 0 ;
  n = cupt -> ngtypes ;
  nvalid = 0 ;
  for (i=0; i<n; i++)  {
   if (cupt->gtypes[i] >= 0) ++nvalid ;  
  }
  return nvalid ;
}  

double malefreq(Indiv *indiv_array, int numindivs) 

/* pop freq of males in sample */
{
   int i ;
   Indiv *indx ; 
   double cmale, cfemale ;

   cmale = 0 ;
   for (i=0; i<numindivs; ++i) { 
    indx = indiv_array + i ;
    if (indx -> gender == 'M') ++cmale ;
   }

   cmale /= (double) numindivs ;

   return cmale ;
}

int isimatch(int a, int b)
{
   if (a < 0)  return YES ;
   if (b < 0)  return YES ;
   if (a==b)   return YES ;
   return NO;
}

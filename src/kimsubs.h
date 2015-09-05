#define MAXPOL 200

#ifndef KIMUTILS
#define TINY 1.0e-8

typedef struct {
  int type ;
  int pop1 ; 
  int pop2 ; 
  int pop3 ; 
  int leaf ;  
  double *mat1 ;
  double *mat2 ;
  int n1 ;
  int n2 ;
  double tau1 ;
  double tau2 ;
} LOGENTRY ; 

typedef struct {
 int *poly ; 
 double *cpoly ;
 int numvars ; 
 int numterms ; 
} POLY ;

typedef struct { 
 double *pars ; 
 int *partype ;
 int *lognums ;
 int npars ;
 int *ssize ;
 POLY *fpoly ;
 double *fprobs ;
 double *dpoly ;
 double *dprobs ;
 LOGENTRY **logs ;
 int numlogs ;
}  PARS  ;

#endif 
#define KIMUTILS

static double ccoeff[MAXPOL+2] ;
static double gegenval[MAXPOL+2] ;
static double jcoeff[MAXPOL+2] ;
static double jacval[MAXPOL+2] ;
static double pco = -1.0, zco = -10, pcmul ;
static double jpco = -1.0, jzco = -10, jpcmul ;

double evalkim(double x, double p, double t, double s2) ;
double evalkim0(double x, double p, double t, double s2) ;
void setccoeff(double p)  ;
void mkgegen(double z)  ;
double gegen(int n, double x) ;
double cfun (int n, double alpha, double x) ;
double jac13(int n,  double x) ;
double jacobi(int n, int a, int b, double x) ;
double norm11(int n) ;
double norm13(int n) ;

double evalpat(double x, double p, double t, double s2) ;
double evalpat0(double x, double p, double t, double s2) ;
void setjcoeff(double p)  ;

void getcc(double *cc, int n, int a, int b) ;
void getcde(double *cc, int n, int a, int b) ;
void getdd(double *cc, int n, int a, int b) ;
int  kindex(int n, int i) ;
double int13(int n, int i) ;
double intx(double y, int kdeg, double t, double s2) ;
double kimint0 (double y, double t, int isupper) ;
double probfix(double y, double t) ;
double getccoeff(int j) ;

// Kimura lifting (moments)
void setccc(double *x, int n) ;
void setcxmat(int n) ;
double *getcxmat(int y) ;
int getcxmatdim() ;
double *getcx(int y) ;
void setggg(double *x, int n) ; 
void getcc(double *cc, int n, int a, int b) ;
void m2t(double *ft, double *fm, int n) ;
void t2m(double *gt, double *gm, int n) ;
void dodrift(double *gt, double *ft, int n, double ttau) ;
void setnn(double *lv, int n) ;
void setlv(double *lv, int n) ;
void mkliftmat(double *mat, int n, double tau) ;
void mkfiltmat(double *mat, int n, int numfilt) ; 
void liftit(double *a, double *b, int n, double tau) ;
int set_lift(int n)  ;
void destroy_lift() ;
double gettau(int a, int c) ;
void printtau() ;
void puttau(int a, int c, double val) ;
int getmaxdeg(int a) ;
int getspecsize() ;
void xspec(int *spec, int speckode) ;
void mkbasespoly(POLY **ppolypt) ;
void setspoly(POLY *pspoly) ;
void print_spoly(POLY *pspoly) ;
void delete_spoly(POLY *pspoly) ;
void load_spoly(POLY *pspoly, int *pp, double *cc, int nv, int nt) ;
int ptrans_spoly(POLY *pspoly, int *ttcount,  int *ad, int *tc) ;  

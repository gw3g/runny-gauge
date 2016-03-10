#include <stdlib.h>
#include <gsl/gsl_math.h>

/*-----------------------------------------------------------------------------------------------*/
/*                                                                "global" over ALL source files */
extern double      Nc;       // colour group SU(Nc)
extern double       g;       // 'the' coupling
extern double    Temp;       // temperature ~ sets relevant units
extern int         Nf;       // active flavours
extern double   beta0;       // qcd beta function [maybe should make as #define??]
extern int        HTL;       // switch for HTL screening [ 0: off, 1: on ]
extern double   kappa;       // coefficient of Debye mass (default = 1.)
extern size_t   calls;       // max number of MC calls
extern int    alf_run;       // running alpha             FLAG
extern double  lambda;       // Lambda / T_c
extern double       J;       // HTL cut
/*-----------------------------------------------------------------------------------------------*/

/* fundamental rep */
#define Cf  (Nc*Nc-1)/(2.*Nc)
#define Tf  (1./2.)
#define df  Nc

/* adjoint rep */
#define Ca  Nc
#define Ta  Nc
#define da  (Nc*Nc-1.)

/*-----------------------------------------------------------------------------------------------*/

typedef enum p_type { 
  B,  // Boson
  F   // Fermion
} p_type;

typedef enum pol { 
  T,  // transverse
  L   // longitudinal
} pol;

double 
   f(double e, p_type X),    // f(e)       for X type
  bf(double e, p_type X);    // \bar{f}(x) for X type

/*-----------------------------------------------------------------------------------------------*/

void qgp(int Nf);           // initialise matrix elements

typedef struct reaction {
  /*
   *  Structure to classify 2-2 QCD reactions. Includes kinematics in terms of
   *  Mendelstam variables and phase space with thermal weights. Matrix elements 
   *  can be tree-level or renormalised.
   */
  int multiplicity;
  p_type particles[4];
  double (*Ms)(double*,double,double);
} reaction;

extern double mD2;                                            // Debye screening mass
extern int nR;                                                // number of reactions
extern reaction *all_R;                                       // pointer to list of reactions

double kernel( double e[4], double s, double t, reaction R);  // thermal weight for a reaction

/*-----------------------------------------------------------------------------------------------*/

/* Monte Carlo integration see "??_prepInt.c" */
extern double lower[5], upper[5];
double integrand(double *, size_t, void *);

/* operators.c */
double 
  xCx( double (*chi)(double,p_type) ),
   xS( double (*chi)(double,p_type) );

/* test function */
struct f_params {double (*chi)(double,p_type);};

/* basis.c */
double eta();

/* x-sec/htl.c */

double  Replace_G(double,double,double,double,double,double);
double *Replace_Q(double,double,double,double,double,double);


#include "core.h"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "cubature.h"

int Nf;         int calls;                 const gsl_rng_type *tt;                     gsl_rng *ws;

/*
 *  Collisional operator    <χ|C|χ>
 *  5d-MC integrator
 *
 *  hep-ph/0302165, (A18)
 */
#include<stdio.h>

/* cubature template */
int c_cub(unsigned dim, const double *x, void *data_,
	   unsigned fdim, double *val)
{
  double *k = (double *)x;
  val[0] = C_integrand( k, (size_t) dim, data_);
  return 0;
}

double xCx( double (*chi)(double,p_type) ) {

  double res, err;
  struct f_params fp = {(*chi)};                                              // chi for inner product
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (5);                       // workspace 5-dim

  if ((double) calls>0) {
    gsl_rng_env_setup (); tt = gsl_rng_default; ws = gsl_rng_alloc (tt);
    gsl_monte_function coll = { C_integrand, 5, &fp };                            // mc function (integrand)
    gsl_monte_vegas_integrate(&coll,lower,upper,5,10000,ws,s,&res,&err);        // perform sampling
    do {
      gsl_monte_vegas_integrate(&coll,lower,upper,5,(size_t)calls,ws,s,&res,&err);
    }
    while (fabs(gsl_monte_vegas_chisq(s)-1.0) > 0.1);
  }
  else {
    double    tol=1e-2;
    int       MaxEvls = -( (int) calls );
    hcubature(1, c_cub, &fp, 5, lower, upper, MaxEvls, 0, tol, ERROR_INDIVIDUAL, &res, &err);
    /*printf("\n -- %g\n\n",res);*/
  }

  if (alf_run) printf("\r  :  %-1.4f  :", Temp/lambda); else printf("\r  :  %03.5f   :", g);

  printf("  %.8f   :    %.5f   :", err/res, gsl_monte_vegas_chisq(s) );
  if ((double) calls>0) {
  gsl_monte_vegas_free(s);                                                    // free memory
  }

  return res;
}

/*
 * preparation of the integrand for the ``source term'' <S|χ>
 *  ( Monte Carlo 1-dim seems to work fine )
 *
 * hep-ph/0010177,  (6.3)
 *
 * 15/08/20 :   working, see "Fsource_CHECK.nb"
 *
 */
double s_integrand(double *args, size_t dim, void *p) {

  /* recover members */
  struct f_params * fp = (struct f_params *)p;
  double (*chi)(double,p_type) = fp->chi;

  /* variable changes */
  double e1   = 1./args[0]-1.,      // args0  ~ e1
         x    = e1/Temp;            // x :=e/T

  double result =                                                             // calculation
    pow(e1,3)*(      da*   f(x,B) * bf(x,B) * chi(x,B)
                + 2.*df*Nf*f(x,F) * bf(x,F) * chi(x,F) );

  result *=  ( 1./pow(args[0],2) )                                            // Jacobian
            *pow(Temp,-1)                                                     // units... T
            *( -1./(M_PI*M_PI) );                                             // prefactors

  return result;
}

double xS( double (*chi)(double, p_type) ) {

  double res, err;
  struct f_params fp = {(*chi)};                                              // chi for inner product
  double sL[1] = {0}, sU[1] = {1};                                            // integration bounds

  gsl_rng_env_setup (); tt = gsl_rng_default; ws = gsl_rng_alloc (tt);
  gsl_monte_function source   = {&s_integrand, 1, &fp};                       // mc function (integrand)
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);                       // workspace in 1 dim
  gsl_monte_vegas_integrate (&source,sL,sU,1,abs(calls),ws,s,&res,&err);           // perform sampling
  gsl_monte_vegas_free (s);                                                   // free memory

  return res;
}

/*
 * Calculate the gain/loss rate
 *
 */

/* cubature template */
int r_cub(unsigned dim, const double *x, void *data_,
	   unsigned fdim, double *val)
{
  double *k = (double *)x;
  val[0] = RATE_integrand( k, (size_t) dim, data_ );
  return 0;
}


double Rate( double e1 ) {

  double res, err;

  double                                                        // boundaries:
  lowerR[4] = {0.+1e-6, 0.+1e-6,  0.,  0.}, 
  upperR[4] = {1.-1e-6, 1.-1e-6,  1.,  M_PI};

  {
    double    tol=1e-2;
    int       MaxEvls = -( (int) calls );
    hcubature(1, r_cub, &e1, 4, lowerR, upperR, MaxEvls, 0, tol, ERROR_INDIVIDUAL, &res, &err);
    /*printf("\n -- %g\n\n",res);*/
  }

  if (alf_run) printf("\r  :  %-1.4f  :", e1/lambda); else printf("\r  :  %03.5f   :", e1/lambda);

  printf("  %.8f   :  ~~   :", err/res  );

  return res;
}


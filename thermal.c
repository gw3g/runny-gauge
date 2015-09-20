#include "core.h"
/* distribution functions */

double fBose( float e ) {
  return 1. / (exp( e ) - 1.);
};

double fFermi( float e ) {
  return 1. / (exp( e ) + 1.);
};

double f(float e, p_type X) {
  switch(X) {
    case F: return fFermi(e);
    case B: return fBose(e);
  }
};

double bf(float e, p_type X) {
  switch(X) {
    case F: return 1.-fFermi(e);
    case B: return 1.+fBose(e);
  }
};

double kernel( double e[4], double s, double t, reaction R) {
  /*
   *  Take a reaction and weight by thermal factors. This is 
   *  the "kernel" of integration.
   *
   *  Output:
   *          Phi({e_i})*|M(s,t,u)|^2
   */
  double phi1234 = 1.;
  for (int i=0; i<2; i++) { 
    phi1234*=f(e[i],R.particles[i]); 
    phi1234*=bf(e[i+2],R.particles[i+2]); 
  }
  return phi1234*(R.Ms(e,s,t))*( (double) R.multiplicity );
};

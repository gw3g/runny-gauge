#include "core.h"
#include <complex.h>
#include <stdio.h>
double J;

/*
 * HTL gluon self-energy
 *
 * hep-ph/0302165, (A11,12)
 *
 * 15/08/15: OK, see "out/plot.nb"
 */
double *Pi(double z, pol X) {
  /*  z := omega/q  */
  double z2=z*z;

  double complex 
    LD  = clog((1.+z)/(1.-z))-I*M_PI, 
    P00 = 1. - 0.5*z*LD,            // Pi_00 =: Pi_L
    Pii = z2 + (1.-z2)*0.5*z*LD;    // Pi_T

  double *P = (double *)malloc( 2*sizeof(double) );

  switch(X) {
    case T:                         // transverse = spatial
      P[0]  = creal( Pii );
      P[1]  = cimag( Pii );         break;
    case L:                         // longitudinal = temporal
      P[0]  = creal( P00 );
      P[1]  = cimag( P00 );         break;
  }

  return P;
}

/*
 * HTL quark self-energy
 *
 * hep-ph/0302165, (A1,2)
 *
 * 15/08/24: OK, see "out/plot.nb"
 */
double *Sig(double z, pol X) {
  /*  z := omega/q  */
  double z2=z*z;

  double complex 
    LD  = clog((1.+z)/(1.-z))-I*M_PI, 
    S0  = 0.5*LD ,                  // Sigma_0
    Si  = -1. + 0.5*z*LD;           // Sigma_i   ( dir. ~ \hat{q} )

  double *S = (double *)malloc( 2*sizeof(double) );

  switch(X) {
    case T:                         // transverse = spatial
      S[0]  = creal( Si );
      S[1]  = cimag( Si );          break;
    case L:                         // longitudinal = temporal
      S[0]  = creal( S0 );
      S[1]  = cimag( S0 );          break;
  }

  return S;                         // return array of [ Re(Sigma), Im(Sigma) ]
}


/*
 * replacement term for bosonic exchange: (s-u)^2/t^2
 *
 * hep-ph/0302165, (A7)
 *
 * 15/08/25: working... NB: don't forget alf_s in numerator!
 * 15/08/29: valgrind suggests memleak here. modifying to mimic fermions  : FIXED
 * 15/09/16: formulae, "doc/gluon_SE.jpg"
 */
double Replace_G(double s, double t, double e3, double e4, double o, double mu2 ) {
  // ---------------------------------------------
  // e3, e4   :   outgoing energies
  // t, o     :   exchanged momenta
  // mu2      :   thermal mass (to multiply HTL)
  // ---------------------------------------------
  double 
    o2=o*o, q2=o2-t, z=o/sqrt(q2), p14p23=(2*e4-o)*(2*e3+o);                   //(e1+e4).(e2+e3)

  /*double m2 = ( Temp*Temp < q2 ) ? 0. : mu2;*/
  double m2 = ( (J*t + Temp*Temp) < 0 ) ? 0. : mu2;

  double complex DeltaL,  DeltaT, repl;

  double *pL = Pi(z, L), *pT = Pi(z, T);

  DeltaL = 1./( q2 + m2*(pL[0] + I*pL[1]) );           free( pL );
  DeltaT = 1./( -t + m2*(pT[0] + I*pT[1]) );           free( pT );

  repl = ( p14p23*DeltaL + (2.*s + t*( 1. + p14p23/q2 ))*DeltaT );

  return repl*conj(repl);
}


/*
 * replacement term for fermionic exchange: s/u, t/u
 *
 * hep-ph/0302165, (A3-5)
 *
 * 15/08/25: working [  free(..)  ]
 * 15/09/16: formulae, "doc/quark_SE.jpg"
 */
double *Replace_Q(double s, double t, double e3, double e4, double o, double mu2) {
  // ---------------------------------------------
  // e3, e4   :   outgoing energies
  // t, o     :   exchanged momenta o = \omega = E4 - E1 ...
  // mu2      :   thermal mass (to multiply HTL)
  // \vec{q}  :   \vec{p4}-\vec{p1} = \vec{p2}-\vec{p3}
  // R        :   P_1 - P_3
  // ---------------------------------------------
  double
    u=-s-t, e1=e4-o, e2=e3+o;

  double                                   // R := exchange mom. (now u-channel)      |
    r0 = -o +e4 - e3,                       // 0-component of R                        |
    r2 = e1*e1+e3*e3+2*e1*e3+u,            // | -\vec{q} + \vec{p4} - \vec{p3} |^2     |
    r  = sqrt(r2),                         // magnitude \vec{r}, direction \hat{r}    |
    z  = r0/r,
    m2 = ( (J*u + Temp*Temp) < 0 ) ? 0. : mu2;
    /*m2 = ( J*r2 > Temp*Temp ) ? 0. : mu2;*/

  double *S0 = Sig(z,L), *Si = Sig(z,T);

  double complex                                              // \cal{R} = R - Sig(R)
    calR0 = r0 - m2*( S0[0] + I*S0[1] )/r,                    // temporal component
    calRi = r  - m2*( Si[0] + I*Si[1] )/r;                    // parallel to \hat{r}

  free(S0);free(Si);

  u*=-1.;
  double complex                                              // four-products  (P_i.\cal{R})
    P1R   = e1*calR0 - (e1*z+0.5*u/r)*calRi +u/2.,            //
    P2R   = e2*calR0 - (e2*z-0.5*u/r)*calRi -u/2.,            //  See "/doc/quark_SE.jpg"
    P4R   = e4*calR0 - (e4*z+0.5*u/r)*calRi +u/2.;            //

  double
    RRs   = -pow(cabs(calR0),2) + pow(cabs(calRi),2),         // \cal{R}.\cal{R}*
    denom = pow(cabs( calR0*calR0 - calRi*calRi), 2),         // denominator...
    num_s = -4.*creal( P1R * conj(P4R) )  - t*RRs,
    num_t = +4.*creal( P1R * conj(P2R) )  - s*RRs;

  //  15/08/25: seems fine
    /*printf("  R.R*    = %.5f \n", RRs);*/
    /*printf("  num_S    = %.5f \n", num_s);*/
    /*printf("  num_T    = %.5f \n", num_t);*/
    /*printf("  |R.R|^2 = %.5f \n", denom);*/

  double *repl = (double *)malloc( 2*sizeof(double) );

  repl[0] = num_s/denom;                                      // replacement: s/u
  repl[1] = num_t/denom;                                      //              t/u

  return repl;
}

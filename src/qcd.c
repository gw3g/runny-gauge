#include "core.h"
#include <stdio.h>

/*-----------------------------------------------------------------------------------------------*/

double                                       Temp,  g,   kappa,  lambda,  beta0,  Nc = 3, mD2, mF2;
int                                                                                   HTL, alf_run;
/*
 * Running coupling, following AP
 *
 *   --> treat Q2 space/time-like differently                             [hep-ph/9411229, 9512336]
 */
double alphas(double Q2) {                                                      // running coupling

  double alf=1., L=log(fabs(Q2)/pow(lambda,2)), alfMax=10.;

    if (Q2>0)                       alf = .5 - atan(L/M_PI)/M_PI;
  else {
    if (Q2<-2.71828*pow(lambda,2))  alf = 1./L;
  }
  alf *= 4*M_PI/beta0;
  return alf_run ? (alf<alfMax?alf:alfMax) : g*g/(4*M_PI);

  /*double alf = (4*M_PI/beta0) / log(fabs(Q2)/pow(lambda,2));*/
  /*return alf_run ? alf : g*g/(4*M_PI);*/
};
/*-----------------------------------------------------------------------------------------------*/

/*
 *  What follows below is a summary of all the 2-2 matrix elements available to a
 *  quark-gluon plasma. Arguments are s,t and IR terms are isolated and screened.
 *  See end of code for definition of qgp(int Nf).
 *                 ___
 *      (1)  -->--/   \-->--  (4)
 *               /     |
 *               | M2  |
 *               |     /
 *      (2)  -->--\___/-->--  (3)
 *
 *  15/08/25: working
 */

double *M_2(double e[4], double s, double t) {          // energies & Mandelstam invariants :

  double  u = -s-t,                                     e1 = e[0],e2 = e[1],e3 = e[2],e4 = e[3],
         s2 = s*s,  u2 = u*u,   t2 = t*t,
         tu = t*u,  us = u*s,   ts = t*s;               // need temp units :
                                                        e1*= Temp;e2*= Temp;e3*= Temp;e4*= Temp;
  double                                                                // parametric dependence
    at = alphas(t),                     at2 = at*at,
    au = alphas(u),                     au2 = au*au,
    as = alphas(s),                     as2 = as*as,
    av = alphas(pow(s*t*u,0.3333)),     av2 = av*av;                    // thermal masses :

  double                                                          mG2 = kappa*mD2*at*pow(Temp,2), 
                                                                  mQ2 = kappa*mF2*au*pow(Temp,2);
  double
    gg1 = 16.*da*Ca*Ca,                         //  = 1152.      >--    pure glue
    qq1 = 8.*(df*df*Cf*Cf)*(1./da),             //  = 16.        \__   quark sector
    qq2 = 16.*(df*Cf)*( Cf-(0.5*Ca) ),          //  = -32./3.    /     group factors
    cc1 = 8.*df*Cf*Cf,                          //  = 128./3.    \__     cross-
    cc2 = 8.*df*Cf*Ca;                          //  = 96.        /      -coupling

  double                                        //  screening replacements, AMY2, p. 13
     rB = Replace_G(s,t,e3,e4,e4-e1,mG2),       //  `B' for boson     [g]
    *rF = Replace_Q(s,t,e3,e4,e4-e1,mQ2),       //  `F' for fermion   [q]
   rt2B = pow(t-mG2,2),                         //  denominator:  t
   ru2F = pow(u-mQ2,2),                         //                u
   rs2B = pow(s,2),       ts2F = pow(s,2),      //  unscreened:   s

  //                                                                    ``dangerous'' terms
                                                                        rB1, rB2, rF1, rF2;
    /*printf("%.5f, %.5f, %.5f, %.5f \n", at,au,mD2,mF2);*/
  // effective mass
  if (!HTL) 
  //    = su/t2                 = (s2+u2)/t2                = s/u                 = t/u
  { rB1 = us/rt2B,          rB2 = (s2+u2)/rt2B,         rF1 = us/ru2F,        rF2 = tu/ru2F;}

  //  HTL replacements
  else 
  //    = su/t2                 = (s2+u2)/t2                = s/u                 = t/u
  { rB1 = (1.-rB)/4.,       rB2 = (1.+rB)/2.,           rF1 = rF[0],          rF2 = rF[1];}

  /*printf("%g\n", rB1);*/

  free(rF);

  double 
    *M2 = (double *)malloc( 7*sizeof(double) );

  //-----------------------------------------------------------------------------------------//
  // Delbrűck :
                                         M2[0] =  gg1*( av2*0. - at2*2.*rB1 - 0.*as2*tu/s2 );
  // (gg <--> gg)
  //-----------------------------------------------------------------------------------------//
  // Møller :
                                                                 M2[1] =  qq1*( at2*rB2 );
  // (q1q2 <--> q1q2)
                                  M2[2] =  qq1*( at2*2.*rB2 ) + qq2*( -au2*s/u -at2*s/t );
  // (q1q1 <--> q1q1)
  //-----------------------------------------------------------------------------------------//
  // Bhabha :
                      M2[3] = qq1*( at2*rB2 +as2*(t2+u2)/s2  ) - qq2*( as2*u/s +at2*u/t );
  // (q1\bar{q1} <--> q1\bar{q1})
                                                          M2[4] =  qq1*( as2*(t2+u2)/s2 );
  // (q1\bar{q1} <--> q2\bar{q2})
  //-----------------------------------------------------------------------------------------//
  // Annihilation :
                                  M2[5] =     cc1*( au2*2.*rF2 ) - cc2*( as2*(t2+u2)/s2 );
  // (q\bar{q} <--> gg)
  //-----------------------------------------------------------------------------------------//
  // Compton :
                                    M2[6] =  -cc1*( as2*u/s + au2*rF1 ) + cc2*( at2*rB2 );
  // (qg <--> qg) 
  //-----------------------------------------------------------------------------------------//
  //  15/08/25:   Checking SIGNS! Ms > 0
  // printf("%.5f,  %.5f,  %.5f,  %.5f,  %.5f,  %.5f,  %.5f\n", 
  //          M2[0], M2[1], M2[2], M2[3], M2[4], M2[5], M2[6]);
  //
  //  15/08/26:   Don't forget the alpha(Q) factors!
  //
  /*if (M2[0]<0) {printf("rB = %.5f, alf = %.5f \n", rB1, at);};*/
  /*printf("%g\n", M2[0]);*/

  for(int i=0; i<7; i++) M2[i] *= pow(4.*M_PI,2);                                       return M2;
}

/*-----------------------------------------------------------------------------------------------*/

int nR;                                                     //                    ... book-keeping
reaction all_Rt[7];                                         //       prepare memory (7 slot array)
reaction *all_R;

double ms(double e[4], double s, double t, int r)   {       // pattern for matrix element function
  double *a = M_2(e,s,t);
  double res = a[r];                                        free(a);                    return res;
};

double ms_d(  double e[4], double s, double t)  {return ms(e,s,t,0);};
double ms_m1( double e[4], double s, double t)  {return ms(e,s,t,1);};
double ms_m2( double e[4], double s, double t)  {return ms(e,s,t,2);};
double ms_b1( double e[4], double s, double t)  {return ms(e,s,t,3);};
double ms_b2( double e[4], double s, double t)  {return ms(e,s,t,4);};
double ms_a(  double e[4], double s, double t)  {return ms(e,s,t,5);};
double ms_c(  double e[4], double s, double t)  {return ms(e,s,t,6);};

/*-----------------------------------------------------------------------------------------------*/

void qgp(int Nf) {
  /*
   *  Prepare table of QGP reactions.
   *
   *  15/06/21: working, but multiplicity factors wrong??   check "~/doc/ms.jpg"
   *  15/09/10: reached consensus on "how to count" with AP...
   */

  nR = 7;                                                       // thermal masses, sans alf & Temp
  mF2 = 4. * M_PI * (Cf/8.)                         ,
  mD2 = 4. * M_PI * (1./3.)*( Ca + Nf*Cf*(df/da) )  ;
  beta0 = (1./3.)*(11.*((double) Nc) - 2.*((double) Nf));       // QCD \beta_0 coeff.

  reaction rr[7] = {
  //{multiplicity,{1,2,3,4}, &ms    }   <---- TEMPLATE

    {1,           {B,B,B,B}, &ms_d  },  //    Delbrűck        (gg <--> gg)
    {8*Nf*(Nf-1), {F,F,F,F}, &ms_m1 },  //    Møller          (q1q2 <--> q1q2)
    {2*Nf,        {F,F,F,F}, &ms_m2 },  //                    (q1q1 <--> q1q1)
    {4*Nf*(Nf-1), {F,F,F,F}, &ms_b2 },  //    Bhabha          (q1\bar{q1} <--> q2\bar{q2})
    {4*Nf,        {F,F,F,F}, &ms_b1 },  //                    (q1\bar{q1} <--> q1\bar{q1})
    {4*Nf,        {F,F,B,B}, &ms_a  },  //    Annihilation    (q\bar{q} <--> gg)
    {8*Nf,        {B,F,B,F}, &ms_c  }   //    Compton         (qg <--> qg)

  };

  for (int i=0;i<nR;i++) { *(all_Rt+i)=*(rr+i); }                                   all_R = all_Rt;
}


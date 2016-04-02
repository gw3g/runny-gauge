#include "core.h"
#include <stdio.h>
reaction *all_R;
int nR;

/* l=2 Legendre function */
double P2 ( double x ) {
  return (0.5)*(3.*x*x - 1.);
}

double Temp, g;

double                                                        // boundaries:
  lower[5] = {0., 0.,  0.,  0.,  0}, 
  upper[5] = {1., 1.,  1.,  1.,  M_PI};

  /*
   *  Returns the "route" structure to be integrated. Here a list of reactions is
   *  inserted and the relevant integrand is constructed, along with dimension and
   *  hypersurface.
   */

double C_integrand_st(double *, size_t , void *);
double C_integrand_qo(double *, size_t , void *);

  /*
   *  Returns the quadratic form to be integrated against the kernel.
   */

double qForm(double c1, double c2, double c3, double c4,  // the chi's
                    double e[4], double s, double t, double u    ) {
  double XXXX;
  /* angles */
  double
    x12 = 1. - s/(2.*e[0]*e[1]),  x13 = 1. + u/(2.*e[0]*e[2]),  x14 = 1. + t/(2.*e[0]*e[3]),
                                  x23 = 1. + t/(2.*e[1]*e[2]),  x24 = 1. + u/(2.*e[1]*e[3]),
                                                                x34 = 1. - s/(2.*e[2]*e[3]);

  /*
   *  (A18): expand (χ1+χ2-χ3-χ4)^2, where χi = c_i X_i
   *  and express X_i X_j in terms of L_ij  using (A19)  
   *
   *  15/09/07 :  See /AP/code/_functionals.c
   */

  /* cos(theta)^2 */
//  return c1*c1  +   2.*c1*c2*x12*x12   -  2.*c1*c3*x13*x13  - 2.*c1*c4*x14*x14
//                +   c2*c2              -  2.*c2*c3*x23*x23  - 2.*c2*c4*x24*x24
//                                       +  c3*c3             + 2.*c3*c4*x34*x34
//                                                            + c4*c4             ;
  /* LegendreP[ 2, x ] */
    XXXX = c1*c1  +   2.*c1*c2*P2(x12)   -  2.*c1*c3*P2(x13)  - 2.*c1*c4*P2(x14)
                  +   c2*c2              -  2.*c2*c3*P2(x23)  - 2.*c2*c4*P2(x24)
                                         +  c3*c3             + 2.*c3*c4*P2(x34)
                                                              + c4*c4             ;
    /*return XXXX/4.;*/
    return 1.;
}

/*-----------------------------------------------------------------------------------------------*/

  /*
   *  Integration variables :         [ the stu-way ]
   *
   *        { s, t, E3, E4, phi }
   */

double C_integrand_st(double *args, size_t dim, void *p) {         // the integrand:

  /* recover members */
  struct f_params * fp = (struct f_params *)p;

  double
    y   = args[0],                                            // args0 = y    ~ s
    z   = args[1],                                            // args1 = z    ~ t
    eps = args[2],                                            // args2 = eps  ~ E3
    xi  = args[3],                                            // args3 = xi   ~ E4
    phi = args[4];                                            // args4 = phi  ~ E1
  double (*chi)(double, p_type) = fp->chi;

  /* variable changes */
  double
    e3  = 1./(eps) - 1. ,
    e4  = 1./(xi)  - 1. ,

  /* kinematic vars */
    s   = y*4.*e3*e4,
    t   = -s*z,
    u   = -s -t,
    b   = -2.*s*( u*e4 + t*e3 ),
    D   = 4.*s*s*t*u*( 4.*e4*e3 - s ),
    e1  = (b+sqrt(D)*cos( phi ))/(2.*s*s),
    e2  = e3 + e4 - e1;

  /* calculation */
  double e[4] = {e1/Temp,e2/Temp,e3/Temp,e4/Temp};

  double result=0.;

  for (int i=0;i<nR;i++) {
    /*R = all_R[i];*/
    /*printf("%.4f \n", (double) all_R[i].multiplicity);*/
    result += kernel(e, s, t, all_R[i])                         // combine reaction kernels
              *qForm( 
                      chi( e[0], all_R[i].particles[0] ),
                      chi( e[1], all_R[i].particles[1] ),
                      chi( e[2], all_R[i].particles[2] ),
                      chi( e[3], all_R[i].particles[3] ),                   e,  s,  t,  u )
             ;
  };

  result *=  ( (4.*e3*e4)*(2.)/pow( eps*xi, 2 ) )               // Jacobian
            *( pow(Temp,2) )                                    // units... T^3
            *( (1./16.)*(1./pow(2.*M_PI, 7)) )/2.     ;         // prefactors

  return result;
};

/*-----------------------------------------------------------------------------------------------*/

  /*
   *  Integration variables :         [ a la AMY ]
   *
   *        { q, omega, E1, E2, phi }
   */

double C_integrand_qo(double *args, size_t dim, void *p) {         // the integrand:

  /* recover members */
  struct f_params * fp = (struct f_params *)p;
  double
    y   = args[0],                                            // args0 = y    ~ q
    z   = args[1],                                            // args1 = z    ~ o := omega
    eps = args[2],                                            // args2 = eps  ~ E1
    xi  = args[3],                                            // args3 = xi   ~ E2
    phi = args[4];                                            // args4 = phi
  double (*chi)(double, p_type) = fp->chi;

  /* variable changes */
  double
    q   = 1./y - 1.,                                          //  + g*T for HARD cutoff
    o   = q*(z),
    e1  = 1./(eps) - 1. + (0.5)*q*(1-z),
    e2  = 1./(xi)  - 1. + (0.5)*q*(1+z),
    e3  = e2 - o,
    e4  = e1 + o;

  /* kinematic vars */
  double
    t   = o*o - q*q,
    s   = (-t/(2*q*q))*(
           ( (e1+e4)*(e2+e3) + q*q ) - cos(phi)*sqrt( (4.*e1*e4 + t)*(4.*e2*e3 + t) )
          ),
    u   = -t-s;

  /* calculation */
  double e[4] = {e1/Temp,e2/Temp,e3/Temp,e4/Temp};

  double result=0.;

  for (int i=0;i<nR;i++) {                                      result += kernel(e, s, t, all_R[i])
                                                                        // combine reaction kernels
              *qForm( chi( e[0], all_R[i].particles[0] ),
                      chi( e[1], all_R[i].particles[1] ),
                      chi( e[2], all_R[i].particles[2] ),
                      chi( e[3], all_R[i].particles[3] ),                          e,  s,  t,  u )
                ;
  };

  result *=  ( q/pow( y*eps*xi, 2 ) )                           // Jacobian
            *( pow(Temp,-1) )                                   // units... T^3
            *( 1./pow(4.*M_PI, 6) )*1.*8./M_PI;                 // prefactors

  /*printf("%g \n", result);*/
  return result;
};

/*-----------------------------------------------------------------------------------------------*/

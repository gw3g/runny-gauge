#include "core.h"
reaction *all_R;
int nR;

/* l=2 Legendre function */
double P2 ( double x ) {
  return (0.5)*(3.*x*x - 1.);
}

double Temp, g;

double                                                        // boundaries:
  lower[5] = {0.,       0.,  0.,  0.,  0}, 
  upper[5] = {1.-1e-6,  1.,  1.,  1.,  M_PI};

  /*
   *  Returns the "route" structure to be integrated. Here a list of reactions is
   *  inserted and the relevant integrand is constructed, along with dimension and
   *  hypersurface.
   *
   *  Integration variables :         [ the stu-way ]
   *
   *        { s, t, E3, E4, phi }
   */

double integrand(double *args, size_t dim, void *p) {         // the integrand:

  /* recover members */
  struct pair * fp = (struct pair *)p;

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

  /* angles */
  double
    x12 = 1. - s/(2.*e1*e2),    x13 = 1. + u/(2.*e1*e3),    x14 = 1. + t/(2.*e1*e4),
                                x23 = 1. + t/(2.*e2*e3),    x24 = 1. + u/(2.*e2*e4),
                                                            x34 = 1. - s/(2.*e3*e4);


  /* calculation */
  double e[4] = {e1/Temp,e2/Temp,e3/Temp,e4/Temp};

  /*
   *  (A18): expand (χ1+χ2-χ3-χ4)^2, where χi = c_i X_i
   *  and express X_i X_j in terms of L_ij  using (A19)  
   *
   *  15/09/07 :  See /AP/code/_functionals.c
   */
  double qForm(double c1, double c2, double c3, double c4) {
  /* cos(theta)^2 */
//  return c1*c1  +   2.*c1*c2*x12*x12   -  2.*c1*c3*x13*x13  - 2.*c1*c4*x14*x14
//                +   c2*c2              -  2.*c2*c3*x23*x23  - 2.*c2*c4*x24*x24
//                                       +  c3*c3             + 2.*c3*c4*x34*x34
//                                                            + c4*c4             ;
  /* LegendreP[2, x] */
    return c1*c1  +   2.*c1*c2*P2(x12)   -  2.*c1*c3*P2(x13)  - 2.*c1*c4*P2(x14)
                  +   c2*c2              -  2.*c2*c3*P2(x23)  - 2.*c2*c4*P2(x24)
                                         +  c3*c3             + 2.*c3*c4*P2(x34)
                                                              + c4*c4             ;
  }

  double result=0.;

  for (int i=0;i<nR;i++) {
    /*R = all_R[i];*/
    /*printf("%.4f \n", (double) all_R[i].multiplicity);*/
    result += kernel(e, s, t, all_R[i])                         // combine reaction kernels
              /*chi( e[3], all_R[i].particles[3]*(*/
                                                  /*chi( e[0], all_R[i].particles[0] )*/
                                                /*+ chi( e[1], all_R[i].particles[1] )*/
                                                /*- chi( e[2], all_R[i].particles[2] )*/
                                                /*- chi( e[3], all_R[i].particles[3] )*/
                  /*)*/
              *qForm( 
                      chi( e[0], all_R[i].particles[0] ),
                      chi( e[1], all_R[i].particles[1] ),
                      chi( e[2], all_R[i].particles[2] ),
                      chi( e[3], all_R[i].particles[3] )  );
  };

  result *=  ( (4.*e3*e4)*(2.)/pow( eps*xi, 2 ) )               // Jacobian
            *( pow(Temp,2) )                                    // units... T^3
            *( (1./16.)*(1./pow(2.*M_PI, 6))*(1./8.) );         // prefactors

  return result;
};


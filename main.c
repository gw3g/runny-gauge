/*
 * Author: greg jackson
 * Date: Sep 17 2015
 * "LO" viscosity, unexpanded
 * 2-2 collisions, HTL screening
 *
 * 15/09/17:  discrepancy with AP for HTL insertions...
 * 15/09/18:  suggests cut-off on small-q for MC
 *
 */

#include "core.h"
#include <stdio.h>

/*-----------------------------------------------------------------------------------------------*/

/* external parameters */     int Nf; double Temp, g, mD2;

                              // ------------------------
                              // (      switches:       )
                              // ------------------------
int             HTL  = 1   ;  // =1 for HTL, =0 for M_eff
double        kappa  = 1.00;  // kappa*mD^2
size_t        calls  = 1e3 ;  // MC calls
int         alf_run  = 1   ;  // =1 for running coupling
double       lambda  = 1.0 ;  // lambda_{QCD}
double            J  = 1.0 ;  // HTL cut

/*-----------------------------------------------------------------------------------------------*/
/*
 *                                              AMY, nll params  [ mu^star, for Nf=1 is a GUESS ] 
 */
double eta1[7] = {27.126, 60.808, 86.47,  106.66, 122.96, 136.38, 147.3 };
double mSt[7]  = {2.765,  2.7,    2.954,  2.957,  2.954,  2.947,  2.940 };
double eta_NLL() { return                                                   // See (4.17) in AMY-II

                          pow(Temp,3.)*eta1[Nf] / (pow(g,4)*log(mSt[Nf]/sqrt(g*g*mD2/(4.*M_PI))));}

void   eval_T(double,double); void   eval_g(double,double); int   points;  // See after main() ...

/*-----------------------------------------------------------------------------------------------*/

int main() {                                        // Main fnc: to explore... T, alpha  dependence

  points = 5;

  for (int nf=0;nf<1;nf++) {                                     // loop over active quark flavours
    Nf = nf; qgp(Nf);
    /*HTL = 0 ; kappa=1.00; eval_g(1e-3,1e1);*/
    /*HTL = 1 ; kappa=1.00; eval_g(1e-3,1e1);*/
    HTL = 1 ; kappa=1.00; eval_T(1.,5.);
    HTL = 0 ; kappa=1.00; eval_T(1.0,5.);
    HTL = 0 ; kappa=0.25; eval_T(1.0,5.);
    /*HTL = 1 ; eval_g(.01,1.);*/
  }

  return 0;
}

/*-----------------------------------------------------------------------------------------------*/

FILE *file; char fname[40];

void eval_T(double Tmin, double Tmax) 
{ alf_run=1; g = 1.; double res1, res2, res3;

       if (!HTL) sprintf(fname, "out/eta(T), M_eff, (kappa=%.2f) Nf=%d.csv", kappa, Nf            );
  else if  (HTL) sprintf(fname, "out/eta(T), HTL, Nf=%d.csv", Nf                                  );

  file = fopen(fname,"w+");

  fprintf(file,   "# GJ, eta w/ only 2->2 processes\n"                                            );
  fprintf(file,   "# Nf = %d, Lambda/Tc=%.3f\n",Nf,lambda                                         );

       if (!HTL) fprintf(file, "# screening: M_eff, kappa = %.3f\n", kappa                        );
  else if  (HTL) fprintf(file, "# screening: htl\n"                                               );

  fprintf(file,   "# 1-fnc basis (Legendre)\n"                                                    );
  fprintf(file,   "# MC samples, %d\n",(int) calls                                                );
  fprintf(file,   "#\n"                                                                           );
  fprintf(file,   "# T,       eta/T^3\n"                                                          );

  printf("\n [ Nf = %d ] \n", Nf );
  printf("  -------------------------------------------------------\n" );
  printf("  : T/lambda :    rel err    :   chisq/dof  :  eta/T^3  :\n" );
  printf("  -------------------------------------------------------\n" );
  for(int i=0; i<points; i++) {
    Temp = Tmin + (Tmax-Tmin)*( ((double) i)/((double) points) );
    printf("  :  %-1.4f  :", Temp/lambda);                       fprintf(file, "%.8f", Temp/lambda);
    J = .5;  res3 = eta()/pow(Temp,3);
    J = 2.;  res2 = eta()/pow(Temp,3);
    J = 1.;  res1 = eta()/pow(Temp,3);
    printf("   %-1.3f   :\n",res1);            fprintf(file, ",%.8f,%.8f,%.8f\n", res1, res2, res3);
  }
  printf("  -------------------------------------------------------\n" );
  fclose(file);
}

void eval_g(double gmin, double gmax) 
{ alf_run=0; Temp = 1.; double res1, res2, res3;

       if (!HTL) sprintf(fname, "out/eta(g), M_eff, (kappa=%.2f) Nf=%d.csv", kappa, Nf            );
  else if  (HTL) sprintf(fname, "out/eta(g), HTL, Nf=%d.csv", Nf                                  );

  file = fopen(fname,"w+");

  fprintf(file,   "# GJ, eta w/ only 2->2 processes\n"                                            );
  fprintf(file,   "# Nf = %d, Lambda/Tc=%.3f\n",Nf,lambda                                         );

       if (!HTL) fprintf(file, "# screening: M_eff, kappa = %.3f\n", kappa                        );
  else if  (HTL) fprintf(file, "# screening: htl\n"                                               );

  fprintf(file,   "# 1-fnc basis (Legendre)\n"                                                    );
  fprintf(file,   "# MC samples, %d\n",(int) calls                                                );
  fprintf(file,   "#\n"                                                                           );
  fprintf(file,   "# g,      eta/T^3\n"                                                           );

  printf("\n [ Nf = %d ] \n", Nf );
  printf("  ---------------------------------------------------------\n" );
  printf("  :    g       :    rel err    :   chisq/dof  :  eta/T^3  :\n" );
  printf("  ---------------------------------------------------------\n" );
  for(int i=0; i<points; i++) {
    g = gmax*pow(10., -(points -1 - i)*( log(gmax/gmin)/log(10.))/((double) points - 1));
    printf("  :  %03.5f   :", g);                                           fprintf(file, "%.8f",g);
    J = .5;  res3 = eta()/pow(Temp,3);
    J = 2.;  res2 = eta()/pow(Temp,3);
    J = 1.;  res1 = eta()/pow(Temp,3);
    printf("  %-1.1e  :\n",res1);            fprintf(file, ",%.8f,%.8f,%.8f\n", res1, res2, res3);
  }
  printf("  ---------------------------------------------------------\n" );
  fclose(file);
}


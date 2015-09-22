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
size_t        calls  = 1e4 ;  // MC calls
int         alf_run  = 1   ;  // =1 for running coupling
double       lambda  = 1.0 ;  // lambda / T_c


/*-----------------------------------------------------------------------------------------------*/
/*
 *                                              AMY, nll params  [ mu^star, for Nf=1 is a GUESS ] 
 */
double eta1[7] = {27.126, 60.808, 86.47,  106.66, 122.96, 136.38, 147.3 };
double mSt[7]  = {2.765,  2.7,    2.954,  2.957,  2.954,  2.947,  2.940 };
double eta_NLL() { return                                                   // See (4.17) in AMY-II

                          pow(Temp,3.)*eta1[Nf] / (pow(g,4)*log(mSt[Nf]/sqrt(g*g*mD2/(4.*M_PI))));}

/*-----------------------------------------------------------------------------------------------*/

int main() {                                  //Main fnc: to explore...        T, alpha  dependence

  int points = 20;                            char fname[40];                           FILE *file;

  for (int nf=0;nf<1;nf++) {                                     // loop over active quark flavours

    Nf=nf;    qgp(Nf);  Temp=1.;  double res;

    /*sprintf(fname, "out/M_eff, (kappa=%.2f) Nf=%d.csv",kappa, Nf);         file = fopen(fname,"w+");*/
    sprintf(fname, "../comparisons/temp.dep/HTL, Nf=%d.csv",kappa, Nf);    file = fopen(fname,"w+");

    fprintf(file,   "# GJ, eta w/ only 2->2 processes\n"                                          );
    fprintf(file,   "# Nf = %d, Lambda/Tc=%.3f\n",Nf,lambda                                       );
    fprintf(file,   "# screening: htl\n"                                                          );
    /*fprintf(file,   "# screening: M_eff, kappa = %.4f\n",kappa                                    );*/
    fprintf(file,   "# 1-fnc basis (Legendre)\n"                                                  );
    fprintf(file,   "# MC samples, %d\n",(int) calls                                              );
    fprintf(file,   "#\n"                                                                         );
    fprintf(file,   "# T,       eta/T^3\n"                                                        );

    printf("\n  :  T/Tc  :    rel err    :   chisq/dof  :  eta/T^3  :\n" );
    printf("  -----------------------------------------------------\n" );
    for(int i=0; i<points; i++) {

      /*g = 10.*pow(10., -(points -1 - i)*4./((double) points - 1));    // evaluate at g = ...*/
      /*printf("%.5f\n", g);*/

      Temp = 1. + 3.*( ((double) i)/((double) points) );
      printf("  :  %-1.2f  :", Temp);
      res = eta()/pow(Temp,3);
      fprintf(file, "%.8f",Temp);  fprintf(file, ",%.8f\n", res );
      printf("   %-1.3f   :\n", res);
    };
    printf("  -----------------------------------------------------\n" );
    fclose(file);
  }

  return 0;
}

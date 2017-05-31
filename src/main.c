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
int           calls  =-1e7 ;  // MC calls {if > 0 : GSL, else hcubature}
int         alf_run  = 0   ;  // =1 for running coupling
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
void   rate_E(double,double,double);

/*-----------------------------------------------------------------------------------------------*/

int main(int argc, char **argv) {                                        // Main fnc: to explore... T, alpha  dependence

  C_integrand = &C_integrand_st;
  points = 60; Temp=1.;Nf=0;

  while (argc--) Nf=(int) atoi(*argv++);
  /*for (int nf=0;nf<1;nf++) {                                     // loop over active quark flavours*/
    /*Nf = nf; */
    qgp(Nf);
    // interaction rate
    /*HTL = 0 ; kappa=1.00; rate_T(.01,9.);*/
    /*HTL = 1 ; kappa=1.00; rate_T(.1,4.,1.0);*/

    /*HTL = 1 ; kappa=1.00; rate_E(0.1,100.,1.0);*/
    /*HTL = 1 ; kappa=1.00; rate_E(1.,8.,0.8);*/
    /*HTL = 1 ; kappa=1.00; rate_E(1.,8.,2.0);*/
    /*HTL = 1 ; kappa=1.00; rate_E(1.,8.,5.0);*/
    /*HTL = 1 ; kappa=1.00; rate_E(1.,20.,2.0);*/
    /*HTL = 1 ; kappa=1.00; rate_E(1.,80.,.8);*/
    /*HTL = 1 ; kappa=1.00; rate_E(1.,200.,2.0);*/
    /*HTL = 1 ; kappa=1.00; rate_E(1.,1000.,10.0);*/

    /*HTL = 1 ; kappa=1.00; rate_T(1.,8.,.8);*/
    /*HTL = 1 ; kappa=1.00; rate_T(.1,4.,3.);*/
    /*HTL = 1 ; kappa=1.00; rate_T(.01,9.);*/
    /*HTL = 0 ; kappa=0.25; rate_T(1e-3,1e0);*/
    /*HTL = 1 ; kappa=1.00; rate_T(1e-3,1e0);*/

    // fixed alpha
    /*HTL = 0 ; kappa=0.25; eval_g(1e-3,1e1);*/
    /*HTL = 0 ; kappa=0.50; eval_g(1e-3,1e1);*/
    /*HTL = 1 ; kappa=1.00; eval_g(1e-3,1e1);*/

    // T-dep
    HTL = 0 ; kappa=0.5;  eval_T(1.0,9.);
    HTL = 1 ; kappa=1.00; eval_T(1.0,9.);

  /*}*/

  return 0;
}

/*-----------------------------------------------------------------------------------------------*/

FILE *file; char fname[90];

void eval_T(double Tmin, double Tmax) 
{ alf_run=1; g = 1.; double res1, res2, res3, res4;

       if (!HTL) sprintf(fname, "out/data/etaT3_kappa%.2f_nf%d_running.dat", kappa, Nf            );
  else if  (HTL) sprintf(fname, "out/data/etaT3_HTL_nf%d_running.dat", Nf                         );

  file = fopen(fname,"w+");

       if (!HTL) fprintf(file, "# screening: M_eff, kappa = %.3f\n", kappa                        );
  else if  (HTL) fprintf(file, "# screening: HTL\n"                                               );

  fprintf(file,   "# MC samples, %d\n",(int) calls                                                );
  fprintf(file,   "#\n"                                                                           );
  fprintf(file,   "# Columns: T/lambda, eta/T^3, {lower, upper, omni}\n"                          );

  printf("\n [ Nf = %d ] ~ ", Nf );
       if (!HTL) printf( " w/ screening via M_eff, kappa = %.3f\n", kappa                        );
  else if  (HTL) printf( " w/ HTL propagators\n"                                                 );
  printf("  --------------------------------------------------------\n" );
  printf("  : T/lambda :    rel err    :   chisq/dof  :  eta/T^3   :\n" );
  printf("  --------------------------------------------------------\n" );
  Temp = 2.;
  for(int i=0; i<points; i++) {
    Temp = Tmin + (Tmax-Tmin)*( ((double) i)/((double) points-1) );
    /*Es = Tmax*pow(10., -(points -1 - i)*( log(Tmax/Tmin)/log(10.))/((double) points - 1));*/
    printf("  :  %-1.4f  :", Temp/lambda);                      fprintf(file, "%.4f", Temp/lambda);
    J = .0;  res4 = eta(); J = .5;  res3 = eta(); J = 2.;  res2 = eta(); J = 1.;  res1 = eta();
    printf("  %1.2e  :\n",res1);fprintf(file, "    %e    %e    %e    %e\n", res1, res2, res3, res4);
  }
  printf("  --------------------------------------------------------\n" );
  fclose(file);
}

void eval_g(double gmin, double gmax) 
{ alf_run=0; Temp = 1.; double res1, res2, res3, res4;

       if (!HTL) sprintf(fname, "out/data/etaT3_kappa%.2f_nf%d_fixed.dat", kappa, Nf              );
  else if  (HTL) sprintf(fname, "out/data/etaT3_HTL_nf%d_fixed.dat", Nf                           );

  file = fopen(fname,"w+");

       if (!HTL) fprintf(file, "# screening: M_eff, kappa = %.3f\n", kappa                        );
  else if  (HTL) fprintf(file, "# screening: HTL\n"                                               );

  fprintf(file,   "# MC samples, %d\n",(int) calls                                                );
  fprintf(file,   "#\n"                                                                           );
  fprintf(file,   "# Columns: g, eta/T^3, {lower, upper, omni} \n"                                );

  printf("\n [ Nf = %d ] ~ ", Nf );
       if (!HTL) printf( " w/ screening via M_eff, kappa = %.3f\n", kappa                         );
  else if  (HTL) printf( " w/ HTL propagators\n"                                                  );
  printf("  ----------------------------------------------------------\n" );
  printf("  :    g       :    rel err    :   chisq/dof  :  eta/T^3   :\n" );
  printf("  ----------------------------------------------------------\n" );
  for(int i=0; i<points; i++) {
    g = gmax*pow(10., -(points -1 - i)*( log(gmax/gmin)/log(10.))/((double) points - 1));
    printf("  :  %03.5f   :", g);                                           fprintf(file, "%.4f",g);
    J = .0;  res4 = eta(); J = .5;  res3 = eta(); J = 2.;  res2 = eta(); J = 1.;  res1 = eta();
    printf("  %1.2e  :\n",res1);fprintf(file, "    %e    %e    %e    %e\n", res1, res2, res3, res4);
  }
  printf("  ----------------------------------------------------------\n" );
  fclose(file);
}

void rate_E(double Emin, double Emax, double TT) 
{ alf_run=1; Temp = TT; double res1, res2, res3, res4; double Ene; 
  //g=sqrt(4.*M_PI*1.);

       if (!HTL) sprintf(fname, "out/data/RT_T%.1f_kappa%.2f_nf%d_running-g.dat", Temp, kappa, Nf  );
  else if  (HTL) sprintf(fname, "out/data/RT_T%.1f_HTL_nf%d_running-g.dat", Temp, Nf               );
       /*if (!HTL) sprintf(fname, "out/data/RT_T%.1f_kappa%.2f_nf%d_fixed-g.dat", Temp, kappa, Nf  );*/
  /*else if  (HTL) sprintf(fname, "out/data/RT_T%.1f_HTL_nf%d_fixed-g.dat", Temp, Nf               );*/

  file = fopen(fname,"w+");

  fprintf(file,   "# Nf = %d, Lambda/Tc=%.3f\n",Nf,lambda                                         );

       if (!HTL) fprintf(file, "# screening: M_eff, kappa = %.3f\n", kappa                        );
  else if  (HTL) fprintf(file, "# screening: htl\n"                                               );

  printf("\n [ Nf = %d ] \n", Nf );
  printf("  ---------------------------------------------------------\n" );
  printf("  :  e/Temp   :    rel err    :             :     R/T     :\n" );
  printf("  ---------------------------------------------------------\n" );
  for(int i=0; i<points; i++) {
    /*g = gmax*pow(10., -(points -1 - i)*( log(gmax/gmin)/log(10.))/((double) points - 1));*/
    /*Ene = Emin + (Emax-Emin)*( ((double) i)/((double) points-1) );*/
    Ene = Emax*pow(10., -(points -1 - i)*( log(Emax/Emin)/log(10.))/((double) points - 1));
    printf("  :  %03.6f   :", Ene/Temp );                        fprintf(file, "%.8f",Ene/Temp);
    J = .0;  res4 =     Rate(Ene)/Temp;
    J = .5;  res3 =     Rate(Ene)/Temp;
    J = 2.;  res2 =     Rate(Ene)/Temp;
    J = 1.;  res1 =     Rate(Ene)/Temp;
    printf("  %-1.3e  :\n",res1);            fprintf(file, "    %e    %e    %e    %e\n", res1, res2, res3, res4);
  }
  printf("  ---------------------------------------------------------\n" );
  fclose(file);
}


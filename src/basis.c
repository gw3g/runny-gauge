#include "core.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define N 1

/*
 * Basis elements, labelled by m = 1,...,N
 * Argument x = p/T is dimensionless!
 *
 * hep-ph/0010177,  (4.9)
 *
 */
double basis(int m, double x) { return pow(x,m+1)/pow(1+x,N-1); }

double func(int n, double x, p_type K) {
  switch (K) {
    case B: return (n>N) ? 0. : basis(n,x);     //  n = 1, 2, ... N           : gluons
    case F: return (n>N) ? basis(n-N,x) : 0.;   //  n = N+1, N+2, ... , 2N    : fermions
  }
}

double C_mn(int m, int n) {
  /* matrix form of collision operator; <m|C|n> */
  double func1(double x, p_type K) {return func(n,x,K)+func(m,x,K);}
  return xCx(&func1);
}

double S_n(int n) {
  /* vector form of "source"; <S|n> */
  double func1(double x, p_type K) {return func(n,x,K);}
  return xS(&func1);
}

double eta() {
  /*
   * Compute the quadratic form; s^t.C^{-1}.s
   * using GSL for vector/matrix classes.
   * result -> proportional to eta
   */

  if (Nf==0) {                                            // Nf = 0, special case
    double ss = S_n(1), cc = C_mn(1,2);
    /*printf("ss=%g, cc=%g \n", ss, cc);*/
    return ss*ss/(15.*cc);
  }
  else {                                                  // Nf > 0, invertible matrix

  double
    *cM = (double *)malloc( (4*N*N) *sizeof(double) ),
    *sV = (double *)malloc( (2*N)   *sizeof(double) );

  for (int i=0; i<2*N; i++) {
      sV[i]=S_n(i+1);
      cM[(i)*2*N+i]=C_mn(i+1,i+1);
  }
  for (int i=0; i<2*N; i++) {
    for (int j=0; j<i; j++) {
      /*
       * Trick to calculate the ``off-diagonal'' elements:
       *
       *      <1|C|2>  =  ( <1+2|C|1+2> - <1|C|1> - <2|C|2> ,
       *
       *  using the linearity of inner prod & symmetry.
       *
       */
      cM[(i)*2*N+j]= 0.5*( C_mn(i+1,j+1) -cM[(i)*2*N+i] -cM[(j)*2*N+j] ) ;
      cM[(j)*2*N+i]= 0.5*( C_mn(i+1,j+1) -cM[(i)*2*N+i] -cM[(j)*2*N+j] ) ;
    }
  }

  gsl_matrix_view cc 
    = gsl_matrix_view_array (cM, 2*N, 2*N);

  /*printf ("M = \n"); gsl_matrix_fprintf (stdout, &cc.matrix, "%g");*/

  gsl_vector_view ss
    = gsl_vector_view_array (sV, 2*N);

  /*printf ("b = \n"); gsl_vector_fprintf (stdout, &ss.vector, "%g");*/

  gsl_vector *x = gsl_vector_alloc (2*N);                 // want to solve for x, where c.x=s

  int s;
  gsl_permutation *p = gsl_permutation_alloc (2*N);       // rand num gen for GSL alg.
  gsl_linalg_LU_decomp (&cc.matrix, p, &s);
  gsl_linalg_LU_solve (&cc.matrix, p, &ss.vector, x);     // determine x, where [C].[x] = [s]

  double qF;
  gsl_blas_ddot(&ss.vector, x, &qF);                      // quadratic form ~ AMY.1  (6.7)
  free(p);free(x);free(cM);free(sV);                      // free mem.
  return qF/15.;
  }
}

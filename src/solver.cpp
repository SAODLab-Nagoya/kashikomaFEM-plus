#include <stdio.h>
#include <math.h>
#include "func_header.h"

/*===== CG method =====*/
void CG_method(double **A, double *x, double *b, int n)
//simultaneous linear equation with n unknowns "Ax = b"
{
  int    i, j;
  int    loop_ctr;     //loop counter
  double alpha, beta;  //coefficients
  double rv[n];        //residual vector
  double rv_old[n];    //rv at previous step
  double pv[n];        //correction DIR_2D vector
  double A_pv[n];      //A * pv
  double b_b;          //b * b
  double err;          //error, rv * rv
  double eps;          //floating-point number called epsilon
  
  #pragma omp parallel for
  for(i = 0; i < n; i++)
  {
    x[i] = 0.0;    //at first, approximate solution x = 0.0
    rv[i] = b[i];  //at first, residual vector r = b - Ax = b
    pv[i] = rv[i];
  }
  b_b = dot_product(b, b, n);
  loop_ctr = 0;
  eps = pow(10, -8);
  
  printf("CG_method started.\n");
  do
  {
    #pragma omp parallel for private(j)
    for(i = 0; i < n; i++)
    {
      A_pv[i] =0.0;
      for (j = 0; j < n; j++)
      {
        A_pv[i] += A[i][j] * pv[j];
      }
    }
    alpha = dot_product(rv, rv, n) / dot_product(pv, A_pv, n);
    #pragma omp parallel for
    for(i = 0; i < n; i++)
    {
      x[i]      += alpha * pv[i];     //update x
      rv_old[i]  = rv[i];             //register current rv
      rv[i]     -= alpha * A_pv[i];   //update rv
    }
    beta = dot_product(rv, rv, n) / dot_product(rv_old, rv_old, n);
    #pragma omp parallel for
    for(i = 0; i < n; i++)
    {
      pv[i] = rv[i] + beta * pv[i];  //update pv
    }
    err = dot_product(rv, rv, n) / b_b;
    printf("loop = %6d, residual = %12.8f\n", ++loop_ctr, err);
  }
  while (err >= eps);
  printf("CG_method finished.\n");
}

/*===== diagonal scaling CG method =====*/
void DSCG_method(double **A, double *x, double *b, int n)
{
  int    i, j;
  int    loop_ctr;      //loop counter
  double alpha, beta;   //coefficients
  double dm[n];         //diagonal component of preprocessing "inverse" matrix
  double rv[n];         //residual vector
  double rv_old[n];     //rv of previous step
  double pv[n];         //correction DIR_2D vector
  double A_pv[n];       //A * pv
  double b_b;           //b * b
  double err;           //error, (rv * rv) / (b * b)
  double eps;           //floating-point number called epsilon
  
  #pragma omp parallel for
  for(i = 0; i < n; i++)
  {
    dm[i] = 1.0 / A[i][i];  //inverce vecter
    x[i]  = 0.0;            //at first, approximate solution x = 0.0
    rv[i] = b[i];           //at first, residual vector r = b - Ax = b
    pv[i] = dm[i] * rv[i];
  }
  b_b = dot_product(b, b, n);
  loop_ctr = 0;
  eps      = pow(10, -8);
  
  printf("DSCG_method started.\n");
  do
  {
    #pragma omp parallel for private(j)
    for(i = 0; i < n; i++)
    {
      A_pv[i] =0.0;
      for (j = 0; j < n; j++)
      {
        A_pv[i] += A[i][j] * pv[j];
      }
    }
    alpha = triple_product(rv, dm, rv, n) / dot_product(pv, A_pv, n);
    #pragma omp parallel for
    for(i = 0; i < n; i++)
    {
      x[i]      += alpha * pv[i];     //update x
      rv_old[i]  = rv[i];             //register current rv
      rv[i]     -= alpha * A_pv[i];   //update rv
    }
    beta = triple_product(rv, dm, rv, n) / triple_product(rv_old, dm, rv_old, n);
    #pragma omp parallel for
    for(i = 0; i < n; i++)
    {
      pv[i] = dm[i] * rv[i] + beta * pv[i];  //update pv
    }
    err = dot_product(rv, rv, n) / b_b;
    printf("loop = %6d, residual = %12.8f\n", ++loop_ctr, err);
  }
  while (err >= eps);
  printf("DSCG_method finished.\n");
}

/*===== calculate inner product of two vectors =====*/
double dot_product(double *vector_a, double *vector_b, int vector_size)
{
  int i;
  double ret = 0.0;  //return value
  
  #pragma omp parallel for reduction(+:ret)
  for(i = 0; i < vector_size; i++)
  {
    ret += vector_a[i] * vector_b[i];
  }
  return ret;
}

/*===== calculate inner product of three vectors =====*/
double triple_product(double *vector_a, double *vector_b, double *vector_c, int vector_size)
{
  int i;
  double ret = 0.0;  //return value
  
  #pragma omp parallel for reduction(+:ret)
  for(i = 0; i < vector_size; i++)
  {
    ret += vector_a[i] * vector_b[i] * vector_c[i];
  }
  return ret;
}

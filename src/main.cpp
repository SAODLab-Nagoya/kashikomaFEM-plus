/*-----------------------------*/
/*  FEM program of c language  */
/*-----------------------------*/

#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "func_header.h"

struct _DATA *data;

/*===== main part =====*/
int main(void)
{
  /*-- for parallel computing --*/
  #ifdef _OPENMP
    omp_set_num_threads(24);
  #endif
  
  /*-- alocate --*/
  data = (DATA*)calloc(1, sizeof(DATA));
  
  /*-- pre_process.c --*/
  input();
  D_matrix();
  Kt_matrix();
  contraction();
  
  /*-- solver.c --*/
  //CG_method(data->Ktc, data->uvc, data->fvc, data->conted_size);
  DSCG_method(data->Ktc, data->uvc, data->fvc, data->conted_size);
  
  /*-- post_process.c --*/
  reversion();
  strain_stress();
  output();
  message();
  
  return 0;
}

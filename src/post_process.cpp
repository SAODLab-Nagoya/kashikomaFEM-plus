#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "func_header.h"

extern struct _DATA *data;

/*===== make uv from uvc =====*/
void reversion()
{
  int i;
  int ctr = 0;  //a counter that shows how many lines to expand
  
  for(i = 0; i < data->dof_total; i++)
  {
    if(data->mark_fixed[i] == 1)  //if dof No.i is fixed
    {
      ctr += 1;
    }
    else
    {
      data->uv[i] = data->uvc[i - ctr];  //reverse operation of contranction
    }
  }
  /*-- free --*/
  free(data->mark_fixed);
  free(data->uvc);
  free(data->fvc);
  free(data->Ktc[0]);
  free(data->Ktc);
}

/*===== calculate strain and stress vector =====*/
void strain_stress()
{
  int i, j;
  int loop_e;  //element loop
  
  /*-- allocate --*/
  data->strain = (double**)calloc(data->nelem, sizeof(double*));
  data->stress = (double**)calloc(data->nelem, sizeof(double*));
  data->strain[0] = (double*)calloc(data->nelem * DIR_2D, sizeof(double));
  data->stress[0] = (double*)calloc(data->nelem * DIR_2D, sizeof(double));
  for(i = 0; i < data->nelem; i++)
  {
    data->strain[i] = data->strain[0] + i * DIR_2D;
    data->stress[i] = data->stress[0] + i * DIR_2D;
  }
  /*-- calculate --*/
  for(loop_e = 0; loop_e < data->nelem; loop_e++)
  {
    /*-- strain = Bm * uv --*/
    for(i = 0; i < DIR_2D; i++)
    {
      for(j = 0; j < DOF_QUAD4; j++)
      {
        data->strain[loop_e][i] += data->ave_Bm[loop_e][i][j] * data->uv[data->connect[loop_e][j/2] * 2 + j % 2];
      }
    }
    /*-- stress = Dm * strain --*/
    for(i = 0; i < DIR_2D; i++)
    {
      for(j = 0; j < DIR_2D; j++)
      {
        data->stress[loop_e][i] += data->Dm[i][j] * data->strain[loop_e][j];
      }
    }
  }
  /*-- free --*/
  free(data->ave_Bm[0][0]);
  free(data->ave_Bm[0]);
  free(data->ave_Bm);
}

/*===== output data to output.vtk =====*/
void output()
{
  int i, j;
  char fname[] = "output_file/output.vtk";
  FILE *output;
  
  /*-- open --*/
  output = fopen(fname, "w");  //write only
  if(output == NULL)
  {  
    printf("%s could not be opened.\n", fname);
    exit(EXIT_FAILURE);
  }
  else
  {
    printf("%s is opened.\n", fname);
  }
  /*--- header ---*/
  fprintf(output, "# vtk DataFile Version 2.0\n");  //first line
  fprintf(output, "output of FEM program\n");       //header
  fprintf(output, "ASCII\n");                       //format
  fprintf(output, "\n");
  fprintf(output, "DATASET UNSTRUCTURED_GRID\n");
  /*--- coordinates of nodes ---*/
  fprintf(output, "POINTS %d float\n", data->nnode);
  for(i = 0; i < data->nnode; i++)
  {
    fprintf(output, "%16.6e%16.6e%16.6e\n", data->xcoord[i], data->ycoord[i], 0.0);
  }
  fprintf(output, "\n");
  /*--- connectivity ---*/
  fprintf(output, "CELLS %d %d\n", data->nelem, (NNODE_QUAD4 + 1) * data->nelem);
  for(i = 0; i < data->nelem; i++)
  {
    fprintf(output, "%8d", NNODE_QUAD4);
    for(j = 0; j < NNODE_QUAD4; j++)
    {
      fprintf(output, "%8d", data->connect[i][j]);
    }
    fprintf(output, "\n");
  }
  fprintf(output, "\n");
  /*--- cell shape ---*/
  fprintf(output, "CELL_TYPES %d\n", data->nelem);
  for(i = 0; i < data->nelem; i++)
  {
    fprintf(output, "%8d\n", 9);  //the cell type number of square is 9
  }
  fprintf(output, "\n");
  /*--- displacement ---*/
  fprintf(output, "POINT_DATA %d\n", data->nnode);
  fprintf(output, "VECTORS displacement float\n");
  for(i = 0; i < data->nnode; i++)
  {
    fprintf(output, "%16.6e%16.6e%16.6e\n", data->uv[2 * i], data->uv[2 * i + 1], 0.0);   //direction x, y, z
  }
  fprintf(output, "\n");
  /*--- strain ---*/
  fprintf(output, "CELL_DATA %d\n", data->nelem);
  fprintf(output, "SCALARS strain_xx float\n");
  fprintf(output, "LOOKUP_TABLE default\n");
  for(i = 0; i < data->nelem; i++)
  {
    fprintf(output, "%16.6e\n", data->strain[i][0]);  //direction xx
  }
  fprintf(output, "\n");
  fprintf(output, "SCALARS strain_yy float\n");
  fprintf(output, "LOOKUP_TABLE default\n");
  for(i = 0; i < data->nelem; i++)
  {
    fprintf(output, "%16.6e\n", data->strain[i][1]);  //direction yy
  }
  fprintf(output, "\n");
  fprintf(output, "SCALARS strain_xy float\n");
  fprintf(output, "LOOKUP_TABLE default\n");
  for(i = 0; i < data->nelem; i++)
  {
    fprintf(output, "%16.6e\n", data->strain[i][2]);  //direction xy
  }
  fprintf(output, "\n");
  /*--- stress ---*/
  fprintf(output, "SCALARS stress_xx float\n");
  fprintf(output, "LOOKUP_TABLE default\n");
  for(i = 0; i < data->nelem; i++)
  {
    fprintf(output, "%16.6e\n", data->stress[i][0]);  //direction xx
  }
  fprintf(output, "\n");
  fprintf(output, "SCALARS stress_yy float\n");
  fprintf(output, "LOOKUP_TABLE default\n");
  for(i = 0; i < data->nelem; i++)
  {
    fprintf(output, "%16.6e\n", data->stress[i][1]);  //direction yy
  }
  fprintf(output, "\n");
  fprintf(output, "SCALARS stress_xy float\n");
  fprintf(output, "LOOKUP_TABLE default\n");
  for(i = 0; i < data->nelem; i++)
  {
    fprintf(output, "%16.6e\n", data->stress[i][2]);  //direction xy
  }
  fprintf(output, "\n");
  fclose(output);
  /*-- free --*/
  free(data->connect[0]);
  free(data->connect);
  free(data->xcoord);
  free(data->ycoord);
  free(data->uv);
  free(data->strain[0]);
  free(data->strain);
  free(data->stress[0]);
  free(data->stress);
}

void message()
{
  /*-- Monospace font recommended --*/
  printf("Output is finished!\n");
  printf("                                                       EE            EE\n");
  printf("      EE           EE             EEEEEEEEEEEL         EE            EF\n");
  printf("      EF    EL     EE                      EF     EEEEEEEEEEEE      EE \n");
  printf("EEEEEEEEEEL  EL    EE                     EF           EE           EF \n");
  printf("     EF   EE  EL   EE                             EEEEEEEEEEEE     EE  \n");
  printf("    EE    EE  IE   EE                                  EE          EF  \n");
  printf("   EF     EE       EE       ET   EL               IEEEEEE              \n");
  printf("  EE     EE        TEL     ET    TEL             EE    EEEEL      EL   \n");
  printf(" EF   TEEF          TEEEEEET      TEEEEEEEEEEE    TEEEEF   TEE    TF   \n");
}

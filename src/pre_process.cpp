#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "func_header.h"

extern struct _DATA *data;

/*===== skip lines of input file =====*/
void dummygets(int loop, FILE *file_pointer)
{
  int  i;
  char buf[256];
  
  for(i = 0; i < loop; i++)
  {
    fgets(buf, sizeof(buf), file_pointer);
  }
}

/*===== input data from input.dat =====*/
void input()
{
  int    a, i, j;
  double tmp;
  char   fname[] = "input_file/input_1.dat";  //name of the input file
  FILE   *input;                              //file pointer
  
  /*-- open --*/
  input = fopen(fname, "r");  //read only
  if(input == NULL)
  {
    printf("%s could not be opened.\n", fname);
    exit(EXIT_FAILURE);  //abnormal termination
  }
  else
  {
    printf("%s is opened.\n", fname);
  }
  /*--- nelem ---*/
  dummygets(2, input);  //skip 2 lines
  fscanf(input, "%d", &data->nelem);
  /*--- nnode, dof_total ---*/
  dummygets(2, input);  //"\n" is remains
  fscanf(input, "%d", &data->nnode);
  data->dof_total = data->nnode * DOF_NODE;
  /*--- connect ---*/
  dummygets(2, input);
  data->connect = (int**)calloc(data->nelem, sizeof(int*));
  data->connect[0] = (int*)calloc(data->nelem * NNODE_QUAD4, sizeof(int));
  for(i = 0; i < data->nelem; i++)
  {
    data->connect[i] = data->connect[0] + i * NNODE_QUAD4;
    fscanf(input, "%d %d %d %d %d", &a, &data->connect[i][0], &data->connect[i][1], &data->connect[i][2], &data->connect[i][3]);
  }
  /*--- xcoord, ycoord ---*/
  dummygets(2, input);
  data->xcoord = (double*)calloc(data->nnode, sizeof(double));
  data->ycoord = (double*)calloc(data->nnode, sizeof(double));
  for(i = 0; i < data->nnode; i++)
  {
    fscanf(input, "%d %lf %lf", &a, &data->xcoord[i], &data->ycoord[i]);
  }
  /*--- uv, fv, conted_size ---*/
  dummygets(3, input);
  fscanf(input, "%d", &data->nfixed);
  data->conted_size = data->dof_total - data->nfixed;  //The size of Ktc is smaller than Kt by nfixed
  dummygets(2, input);
  data->uv = (double*)calloc(data->dof_total, sizeof(double));
  data->mark_fixed = (int*)calloc(data->dof_total, sizeof(int));
  for(i = 0; i < data->nfixed; i++)
  {
    fscanf(input, "%d %lf", &j, &tmp);
    data->uv[j] = tmp;
    data->mark_fixed[j] = 1;  //mark with 1
  }
  dummygets(3, input);
  fscanf(input, "%d", &a);
  dummygets(2, input);
  data->fv = (double*)calloc(data->dof_total, sizeof(double));
  for(i = 0; i < a; i++)
  {
    fscanf(input, "%d %lf", &j, &tmp);
    data->fv[j] = tmp;
  }
  /*--- Young, Poiss ---*/
  dummygets(3, input);
  fscanf(input, "%lf", &data->Young);
  dummygets(2, input);
  fscanf(input, "%lf", &data->Poiss);
  fclose(input);  //close the input file
}

/*===== calculate D matrix =====*/
void D_matrix()
{
  if(PLANE_STRAIN_STATE == 0)
  {
    printf("Calculate in plane strain state.\n");
    double coef = data->Young / ((1 + data->Poiss) * (1 - data->Poiss));   //coefficient
    data->Dm[0][0] = coef;
    data->Dm[0][1] = coef * data->Poiss;
    data->Dm[0][2] = 0;
    data->Dm[1][0] = data->Dm[0][1];
    data->Dm[1][1] = coef;
    data->Dm[1][2] = 0;
    data->Dm[2][0] = data->Dm[0][2];
    data->Dm[2][1] = data->Dm[1][2];
    data->Dm[2][2] = coef * (1 - data->Poiss) / 2;
  }
  else
  {
    printf("Calculate in plane stress state.\n");
    double coef = (data->Young * (1 - data->Poiss)) / ((1 + data->Poiss) * (1 - 2 * data->Poiss));  //coefficient
    data->Dm[0][0] = coef;
    data->Dm[0][1] = coef * data->Poiss / (1 - data->Poiss);
    data->Dm[0][2] = 0;
    data->Dm[1][0] = data->Dm[0][1];
    data->Dm[1][1] = coef;
    data->Dm[1][2] = 0;
    data->Dm[2][0] = data->Dm[0][2];
    data->Dm[2][1] = data->Dm[1][2];
    data->Dm[2][2] = coef * (1 - 2 * data->Poiss) / (2 * (1 - data->Poiss));
  }
}

/*===== calculate Kt matrix =====*/
void Kt_matrix()
{
  int    i, j, k;
  int    loop_e, loop_n, loop_g;                        //element, node, gauss point loop
  int    row_Kt, column_Kt;                             //index number of Kt matrix
  double x[NNODE_QUAD4], y[NNODE_QUAD4];                //coordinates of nodes
  double xi, et, dxdxi, dydxi, dxdet, dydet;            //xi, eta, derivatives
  double dNdx[NGAUSS], dNdy[NGAUSS];                    //derivative of the shape function N
  double gauss_xicoord[NGAUSS], gauss_etcoord[NGAUSS];  //coordinates of gauss points
  double weight_xi[NGAUSS], weight_et[NGAUSS];          //weights of gauss points
  double detJ;                                          //Jacobian
  double Bm[DIR_2D][DOF_QUAD4];                         //B matrix
  double BmtDm[DIR_2D][DOF_QUAD4];                      //Bm^T * Dm
  double BmtDmBm[DOF_QUAD4][DOF_QUAD4];                 //Bm^T * Dm * Bm
  double Ke[DOF_QUAD4][DOF_QUAD4];                      //element stiffness matrix
  
  gauss_xicoord[0] = - pow(3, -0.5);
  gauss_xicoord[1] =   pow(3, -0.5);
  gauss_xicoord[2] = - pow(3, -0.5);
  gauss_xicoord[3] =   pow(3, -0.5);
  gauss_etcoord[0] = - pow(3, -0.5);
  gauss_etcoord[1] = - pow(3, -0.5);
  gauss_etcoord[2] =   pow(3, -0.5);
  gauss_etcoord[3] =   pow(3, -0.5);
  weight_xi[0] = 1.0;
  weight_xi[1] = 1.0;
  weight_xi[2] = 1.0;
  weight_xi[3] = 1.0;
  weight_et[0] = 1.0;
  weight_et[1] = 1.0;
  weight_et[2] = 1.0;
  weight_et[3] = 1.0;
  
  /*-- allocate ave_Bm --*/
  data->ave_Bm = (double***)calloc(data->nelem, sizeof(double**));
  data->ave_Bm[0] = (double**)calloc(data->nelem * DIR_2D, sizeof(double*));
  data->ave_Bm[0][0] = (double*)calloc(data->nelem * DIR_2D * DOF_QUAD4, sizeof(double));
  for(i = 0; i < data->nelem; i++)
  {
    data->ave_Bm[i] = data->ave_Bm[0] + i * DIR_2D;
    for(j = 0; j < DIR_2D; j++)
    {
      data->ave_Bm[i][j] = data->ave_Bm[0][0] + (i * DIR_2D + j) * DOF_QUAD4;
    }
  }
  
  /*-- allocate Km --*/
  data->Kt = (double**)calloc(data->dof_total, sizeof(double*));
  data->Kt[0] = (double*)calloc(data->dof_total * data->dof_total, sizeof(double));
  for(i = 0; i < data->dof_total; i++)
  {
    data->Kt[i] = data->Kt[0] + i * data->dof_total;
  }
  
  /*-- calculate B_matrix --*/
  for(loop_e = 0; loop_e < data->nelem; loop_e++)  //start element loop
  {
    for(i = 0; i < DOF_QUAD4; i++)
    {
      for(j = 0; j < DOF_QUAD4; j++)
      {
        Ke[i][j] = 0.0;  //Ke is determined for each element
      }
    }
    for(loop_n = 0; loop_n < NNODE_QUAD4; loop_n++)
    {
      x[loop_n] = data->xcoord[data->connect[loop_e][loop_n]];
      y[loop_n] = data->ycoord[data->connect[loop_e][loop_n]];
    }
    for(loop_g = 0; loop_g < NGAUSS; loop_g++)  //start gauss point loop
    {
      xi = gauss_xicoord[loop_g];
      et = gauss_etcoord[loop_g];
      dxdxi = (- (1-et) * x[0] + (1-et) * x[1] + (1+et) * x[2] - (1+et) * x[3]) / 4;
      dydxi = (- (1-et) * y[0] + (1-et) * y[1] + (1+et) * y[2] - (1+et) * y[3]) / 4;
      dxdet = (- (1-xi) * x[0] - (1+xi) * x[1] + (1+xi) * x[2] + (1-xi) * x[3]) / 4;
      dydet = (- (1-xi) * y[0] - (1+xi) * y[1] + (1+xi) * y[2] + (1-xi) * y[3]) / 4;
      detJ = dxdxi * dydet - dxdet * dydxi;
      dNdx[0] = (- (1-et) * dydet + (1-xi) * dydxi) / (4 * detJ);
      dNdx[1] = (  (1-et) * dydet + (1+xi) * dydxi) / (4 * detJ);
      dNdx[2] = (  (1+et) * dydet - (1+xi) * dydxi) / (4 * detJ);
      dNdx[3] = (- (1+et) * dydet - (1-xi) * dydxi) / (4 * detJ);
      dNdy[0] = (  (1-et) * dxdet - (1-xi) * dxdxi) / (4 * detJ);
      dNdy[1] = (- (1-et) * dxdet - (1+xi) * dxdxi) / (4 * detJ);
      dNdy[2] = (- (1+et) * dxdet + (1+xi) * dxdxi) / (4 * detJ);
      dNdy[3] = (  (1+et) * dxdet + (1-xi) * dxdxi) / (4 * detJ);
      Bm[0][0] = dNdx[0];
      Bm[0][1] = 0.0;
      Bm[0][2] = dNdx[1];
      Bm[0][3] = 0.0;
      Bm[0][4] = dNdx[2];
      Bm[0][5] = 0.0;
      Bm[0][6] = dNdx[3];
      Bm[0][7] = 0.0;
      Bm[1][0] = 0.0;
      Bm[1][1] = dNdy[0];
      Bm[1][2] = 0.0;
      Bm[1][3] = dNdy[1];
      Bm[1][4] = 0.0;
      Bm[1][5] = dNdy[2];
      Bm[1][6] = 0.0;
      Bm[1][7] = dNdy[3];
      Bm[2][0] = dNdy[0];
      Bm[2][1] = dNdx[0];
      Bm[2][2] = dNdy[1];
      Bm[2][3] = dNdx[1];
      Bm[2][4] = dNdy[2];
      Bm[2][5] = dNdx[2];
      Bm[2][6] = dNdy[3];
      Bm[2][7] = dNdx[3];
      for(i = 0; i < DIR_2D; i++)
      {
        for(j = 0; j < DOF_QUAD4; j++)
        {
          data->ave_Bm[loop_e][i][j] += Bm[i][j] / NGAUSS;  //average Bm of 4 gauss points in each element
        }
      }
      /*-- calculate Ke_matrix --*/
      for(i = 0; i < DOF_QUAD4; i++)
      {
        for(j = 0; j < DIR_2D; j++)
        {
          BmtDm[j][i] = 0.0;  //initialize
          for(k = 0; k < DIR_2D; k++)
          {
            BmtDm[j][i] += Bm[k][i] * data->Dm[k][j];
          }
        }
        for(j = 0; j < DOF_QUAD4; j++)
        {
          BmtDmBm[i][j] = 0.0;  //initialize
          for(k = 0; k < DIR_2D; k++)
          {
            BmtDmBm[i][j] += BmtDm[k][i] * Bm[k][j];
          }
          Ke[i][j] += BmtDmBm[i][j] * detJ * weight_xi[loop_g] * weight_et[loop_g] * THICK;
        }
      }
    }  //end of gauss point loop
    /*-- calculate Kt_matrix --*/
    for(i = 0; i < DOF_QUAD4; i++)
    {
      row_Kt = data->connect[loop_e][i/2] * 2 + i % 2;  //ex. node No.8 and direction y --> row_Kt = 17
      for(j = 0; j < DOF_QUAD4; j++)
      {
        column_Kt = data->connect[loop_e][j/2] * 2 + j % 2;
        data->Kt[row_Kt][column_Kt] += Ke[i][j];
      }
    }
  }  //end of element loop
}

/*===== make contracted matrix and vector =====*/
void contraction()
{
  int i, r, c;
  int ctr_row, ctr_column;  //a counter that shows how many lines to shorten
  
  /*-- allocate --*/
  data->uvc = (double*)calloc(data->conted_size, sizeof(double));
  data->fvc = (double*)calloc(data->conted_size, sizeof(double));
  data->Ktc = (double**)calloc(data->conted_size, sizeof(double*));
  data->Ktc[0] = (double*)calloc(data->conted_size * data->conted_size, sizeof(double));
  for(i = 0; i < data->conted_size; i++)
  {
    data->Ktc[i] = data->Ktc[0] + i * data->conted_size;
  }
  /*-- contraction --*/
  ctr_row = 0;
  for(r = 0; r < data->dof_total; r++)  //row loop
  {
    ctr_column = 0;
    if(data->mark_fixed[r] == 1)  //if dof No.r is fixed
    {
      ctr_row++;
    }
    else
    {
      data->uvc[r - ctr_row] = data->uv[r];
      data->fvc[r - ctr_row] = data->fv[r];
      for(c = 0; c < data->dof_total; c++)  //column loop
      {
        if(data->mark_fixed[c] == 1)
        {
          ctr_column++;
        }
        else
        {
          data->Ktc[r - ctr_row][c - ctr_column] = data->Kt[r][c];
        }
      }
    }
  }
  /*-- free --*/
  free(data->fv);
  free(data->Kt[0]);
  free(data->Kt);
}

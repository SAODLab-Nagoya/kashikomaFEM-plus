#ifndef FUNC_HEADER
#define FUNC_HEADER

/*-- macro --*/
#define DOF_NODE 2            //dof in a node
#define NNODE_QUAD4 4         //number of nodes in a quad4 element
#define DOF_QUAD4 8           //dof in a quad4 element
#define DIR_2D 3              //xx, yy, xz dirctions, in 2D structure
#define THICK 1.0             //thickness
#define NGAUSS 4              //number of gauss points in a element
#define PLANE_STRAIN_STATE 0  //= 0 plane stress state, != 0 plane strain state

/*-- structure --*/
typedef struct _DATA
{
  int    nelem;               //number of total elements
  int    nnode;               //number of total nodes
  int    dof_total;           //number of total dof
  int    **connect;           //node No. = connect[element No.][1, 2, 3, 4]
  double *xcoord, *ycoord;    //coordinate = xcoord[node No.]
  int    nfixed;              //number of fixed dof
  int    *mark_fixed;         //mark_fixed[dof No.] = 1 fixed, = 0 not fixed 
  double *uv, *fv;            //displacement vector, force vector
  double *uvc, *fvc;          //contracted vector
  int    conted_size;         //size of contracted vectors or matrix
  double Young, Poiss;        //Young's modulus, Poisson's ratio
  double Dm[DIR_2D][DIR_2D];  //D or C matrix
  double ***ave_Bm;           //average Bm of 4 gauss points in each element, used in strain calculation
  double **Kt;                //total stiffness matrix
  double **Ktc;               //contracted Kt matrix
  double **strain;            //strain in each element, strain[elem No.][xx, yy, xy]
  double **stress;            //stress in each element, stress[elem No.][xx, yy, xy]
}
DATA;

/*-- pre_process.c --*/
void dummygets(int loop, FILE *file_pointer);
void input();
void D_matrix();
void Kt_matrix();
void contraction();

/*-- post_process.c --*/
void reversion();
void strain_stress();
void output();
void message();

/*-- solver.c --*/
void   CG_method(double **A, double *x, double *b, int n);
void   DSCG_method(double **A, double *x, double *b, int n);
double dot_product(double *vector_a, double *vector_b, int vector_size);
double triple_product(double *vector_a, double *vector_b, double *vector_c, int vector_size);

#endif

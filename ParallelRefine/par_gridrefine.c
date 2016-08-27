#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <mpi.h>

#include "grid_utils.h"
#include "interface_comms.h"
#include "compute_timestep.h"

int rank, size, cRank, cRank2;
MPI_Comm MPI_COMM_CART;
MPI_Comm MPI_COMM_CART2;
int P, Q;
int crds [2];
int rcrds[2];
int *coarse_ranks;
int *refine_ranks;
MPI_Datatype MPI_CLM;
MPI_Datatype MPI_CLM_FINE;


//updates matrices
void update(double*** mat1_ptr, double*** mat2_ptr)
{
  double **tmp = *mat1_ptr;
  *mat1_ptr = *mat2_ptr;
  *mat2_ptr = tmp;
}

int main(int argc, char *argv[])
{
  int i, j, c;
  double **mat1, **mat2;

  //size of grid
  int nrows = 10, ncols = 10;

  //2x2 for mesh in a 4x4 array
  double dx = 0, dy= 0, dt, time;
  int iter, max_iter = 100;

  //MPI Init
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //dimensions of ranks on coarse/fine grid
  Q = 2;
  P = 2;

  while ((c = getopt (argc, argv, "r:i:")) != -1)
    {
      switch(c)
  	{
  	case 'r':
  	  nrows = ncols = atoi(optarg); break;
  	case 'i':
  	  max_iter = atoi(optarg); break;
  	default:
  	  fprintf(stderr, "Invalid option\n");
  	  return -1;
  	}
    }


  // find out where all the raks are going to go - needed for interface comms

  int *initial_coarse_ranks = calloc(P*Q, sizeof(int));	  
  coarse_ranks = calloc(P*Q, sizeof(int));
  refine_ranks = calloc(P*Q, sizeof(int));

  int coarse_index=0, refine_index=0;

  for(i=0; i<Q; i++)
    {
    for(j=0; j<P; j++)
      {
	initial_coarse_ranks[i*P + j] = i*P + j;

  	if(i >= Q/2 && j < P/2) {
  	  refine_ranks[refine_index++] = i*P + j;
	} else {
  	  coarse_ranks[coarse_index++] = i*P + j;
	}
      }
    }
  
  for(i=P*Q; i<size; i++)
    {
      refine_ranks[refine_index++] = i;
    }

  if(!rank) {
    printf("Coarse:\n");
    for(i=0; i<0.75*P*Q; i++)
      printf("%d\n", coarse_ranks[i]);

    printf("Refine:\n");
    for(i=0; i<P*Q; i++)
      printf("%d\n", refine_ranks[i]);
  }

  //parallel rows & cols coarse grid
  int par_rows = 2 + nrows/Q;
  int par_cols = 2 + ncols/P;

  //size of fine grid
  int szofmesh = .5*nrows;
  int refrows = szofmesh * 2 - 1;
  int refcols = szofmesh * 2 - 1;

  //parallel refined rows & cols
  int par_ref_rows = 2 + refrows/Q;
  int par_ref_cols = 2 + refcols/P;

  MPI_Group coarse_group;
  MPI_Group refine_group;
  MPI_Group world_group;

  MPI_Comm coarse_comm;
  MPI_Comm refine_comm;

  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  MPI_Group_incl(world_group, 4, initial_coarse_ranks, &coarse_group);
  MPI_Comm_create(MPI_COMM_WORLD, coarse_group, &coarse_comm);

  MPI_Group_incl(world_group, 1, refine_ranks, &refine_group);
  MPI_Comm_create(MPI_COMM_WORLD, refine_group, &refine_comm);

  int in_refine = !(rank < P*Q);

  MPI_Type_vector(par_rows - 2, 1, par_cols, MPI_DOUBLE, &MPI_CLM);
  MPI_Type_commit(&MPI_CLM);
  MPI_Type_vector(par_ref_rows - 2, 1, par_ref_cols, MPI_DOUBLE, &MPI_CLM_FINE);
  MPI_Type_commit(&MPI_CLM_FINE);

  if(!in_refine) {
    int prds[] = {0,0};
    int dims[] = {Q,P};
    MPI_Cart_create(coarse_comm, 2, dims, prds, 0, &MPI_COMM_CART);
    MPI_Comm_rank(MPI_COMM_CART, &cRank);
    MPI_Cart_coords(MPI_COMM_CART, cRank, 2, crds);
  } 

  for(i=0; i<Q; i++)
    {
    for(j=0; j<P; j++)
      {
  	if(i >= Q/2 && j < P/2) 
        {
          in_refine = i*P + j;
	} 
        else 
        {
  	  coarse_ranks[coarse_index++] = i*P + j;
	}
      }
    }
 
  if (in_refine)
  {
    //Create refine cart
    int rprds[] = {0,0};
    int rdims[] = {Q,P};
    MPI_Cart_create(refine_comm, 2, rdims, rprds, 0, &MPI_COMM_CART2);
    MPI_Comm_rank(MPI_COMM_CART2, &cRank2);
    MPI_Cart_coords(MPI_COMM_CART2, cRank2, 2, rcrds);
  }


  for(i = Q * P; i < size; i++)
    in_refine; 

  if(!in_refine) {
    rcrds[0]=-1;
    rcrds[1]=-1;
  } else {
    crds[0]=-1;
    crds[1]=-1;
  }
 
  if(in_refine) {
    for(i=0; i<Q; i++) {
      for(j=0; j<P; j++) {
	if(rcrds[0] == j && rcrds[1] == i) {
	  printf("%d\n", rank);
	}
      }
    }
  }
  
  //width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  dt = .001;

  //sets up matrices for refinement and not in refinement
  if(!in_refine)
    {
      mat1 = calloc(par_rows, sizeof(double *));
      mat2 = calloc(par_rows, sizeof(double *));

      double *mat1_1d = calloc(par_rows*par_cols, sizeof(double));
      for (i = 0; i < par_rows; i++)
	{
	  mat1[i] = &(mat1_1d[i * par_cols]);
	}
	
      double *mat2_1d = calloc(par_rows*par_cols, sizeof(double));
      for (i = 0; i < par_rows; i++)
	{
	  mat2[i] = &(mat2_1d[i * par_cols]);
	}
      initGrid(mat1, par_rows, par_cols);
      initGrid(mat2, par_rows, par_cols);
    }
  else
    {
      mat1 = calloc(par_ref_rows, sizeof(double *));
      mat2 = calloc(par_ref_rows, sizeof(double *));

      double *mat1_1d = calloc(par_ref_rows*par_ref_cols, sizeof(double));
      for (i = 0; i < par_ref_rows; i++)
	{
	  mat1[i] = &(mat1_1d[i * par_ref_cols]);
	}
	
      double *mat2_1d = calloc(par_ref_rows*par_ref_cols, sizeof(double));
      for (i = 0; i < par_ref_rows; i++)
	{
	  mat2[i] = &(mat2_1d[i * par_ref_cols]);
	}
      initFineGrid(mat1, par_ref_rows, par_ref_cols);
      initFineGrid(mat2, par_ref_rows, par_ref_cols);
    }

  //buffers to send data
  double *bfr_right_refine_pts, *bfr_top_refine_pts;
  bfr_right_refine_pts = calloc(par_ref_rows - 2, sizeof(double));
  bfr_top_refine_pts = calloc(par_ref_cols - 2, sizeof(double));

  double *bfr_right_interface_pts, *bfr_top_interface_pts;
  bfr_right_interface_pts = calloc(par_rows - 2, sizeof(double));
  bfr_top_interface_pts = calloc(par_cols - 2, sizeof(double));


  time = 0;

  //weights
  double wx = 0, wy = 0;
  //diffusitivity coefficient
  double k = .02;

  //printf("dx = %lf, dy = %lf, dt = %lf\n", dx, dy, dt);
  //weight in x direction
  wx = k * dt/(dx*dx);
  //weight in y direction
  wy = k * dt/(dy*dy);

  //stability
  // wx + wy must be less than .5
  if ( wx + wy > .5)
    {
      printf("Does not reach requirements for stability\n");
      printf("wx = %lf, wy = %lf\n", wx, wy);
      exit(1);
    }

  double wxy = k * dt/(dx*dx);
  int refinement = 2;
  double sxy = (double)dx/refinement;
  
  for(iter = 0; iter < max_iter; iter++)
    {
      time = time + dt;

      if(in_refine || rank == 3)
        sendRightPoints(mat1, bfr_right_interface_pts, par_rows, par_cols);

      if(in_refine || rank == 0)
        sendTopPoints(mat1, bfr_top_interface_pts, par_rows, par_ref_cols);

      if(in_refine)
      {
	 injectFinePoints(bfr_top_interface_pts, bfr_right_interface_pts, mat1, par_ref_rows, par_cols, par_ref_cols);
	 computeInterfaceTopTimestep(bfr_top_refine_pts, mat1, mat2, par_rows, par_ref_cols, wxy, dx);
	 computeInterfaceRightTimestep(bfr_right_refine_pts, mat1, mat2, par_ref_rows, ncols, wxy, dx, refcols);
	 computeCornerTimestep(bfr_right_refine_pts, mat1, mat2, wxy, par_ref_cols);
      }

      if(in_refine || rank == 3)
	sendRightInterfacePoints(mat1, bfr_right_refine_pts, par_ref_rows, par_ref_cols);

      if(in_refine || rank == 0)
	sendTopInterfacePoints(mat1, bfr_top_refine_pts, par_ref_rows, par_ref_cols);

      if (rank == 3 || rank == 0)
        injectCoarsePoints(bfr_top_refine_pts, bfr_right_refine_pts, mat1, par_rows, par_cols);


      if(!in_refine){

	computeTimestep(bfr_top_refine_pts, bfr_right_refine_pts, mat1, mat2, par_rows, par_cols, wx, wy);
      }
      else
      {
	    computeFineTimestep(mat1, mat2, refrows, refcols, wxy, sxy);
      }

	update(&mat1, &mat2);

    }

  if(!in_refine)
    printGrid(mat1, par_rows, par_cols);
  else
  {
    printGrid(mat1, par_ref_rows, par_ref_cols);
  }

  MPI_Finalize();

  return 0;
}


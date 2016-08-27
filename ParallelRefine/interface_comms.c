#include "interface_comms.h"
#include <mpi.h>

extern int P, Q, rank, size;
extern int crds[2];
extern int rcrds[2];
extern int *coarse_ranks;
extern int *refine_ranks;
extern MPI_Comm MPI_COMM_CART, MPI_COMM_CART2;

//injects buffers into coarse grid
void injectCoarsePoints(double* toparray, double* rightarray, double** mat, int par_rows, int par_cols)
{
  int i, j = 0;
  //top
  if (rank == 0){
  for (i = 1; i < par_cols - 1; i++)
  {
    mat[par_rows - 1][i] = toparray[j];
    j+=2;
  }
  }
  //right
  if (rank == 3){
  for(i = 1; i < par_rows - 1; i++)
    {
      mat[i][0] = rightarray[j];
      j+=2;
    }
  }
}

//injects buffers into refined grid
void injectFinePoints(double* topbuffer, double* rightbuffer, double** mat, int par_ref_rows, int par_cols, int par_ref_cols)
{
  int i, j = 0, n = 0;

  //top
  for (i = 1; i < par_ref_cols - 1; i+=2)
    {
      mat[0][i] = topbuffer[j];
      j++;
    }
  //right
  for(i = 1; i < par_ref_rows - 1; i+=2)
    {
      mat[i][par_ref_cols - 1] = rightbuffer[n];
      n++;
    }

  //averages points on fine grid
  for(i=2; i<par_ref_cols-2; i+=2)
    {
      mat[1][i] = 0.5*(mat[1][i-1] + mat[1][i+1]);
    }

  for(i=2; i<par_ref_rows-2; i+=2)
    {
      mat[i][par_ref_cols-2] = 0.5*(mat[i-1][par_ref_cols-2] + mat[i+1][par_ref_cols-2]);
    }
}

//sends points from coarse grid to right of refined grid
void sendRightPoints(double** mat, double* buffer, int par_rows, int par_cols)
{
  MPI_Status stat[2]; 
  MPI_Datatype MPI_CLM_BFR;
  MPI_Type_vector(0.5*(par_rows - 2), 1, par_cols, MPI_DOUBLE, &MPI_CLM_BFR);
  MPI_Type_commit(&MPI_CLM_BFR);

  if(crds[0] == P/2 && crds[1] >= Q/2) {

    int first_target  = refine_ranks[(   Q/2) *P*(crds[1]-(Q/2)+1)];
    int second_target = refine_ranks[(1+(Q/2))*P*(crds[1]-(Q/2)+1)];

    MPI_Send(&mat[1][1], 1, MPI_CLM_BFR, first_target, 0, MPI_COMM_WORLD);

    MPI_Send(&mat[1][par_rows/2], 1, MPI_CLM_BFR, second_target, 0, MPI_COMM_WORLD);
  }

  if(rcrds[0] == P-1) {
    MPI_Recv(buffer, 0.5*(par_rows-2), MPI_DOUBLE, coarse_ranks[(Q/2)*P + (P/2)*(rcrds[0]/2)], 0, MPI_COMM_WORLD, &stat[1]);
  }

}

//sends points from coarse grid to top of refined grid
void sendTopPoints(double** mat, double* buffer, int par_rows, int par_ref_cols)
{
  MPI_Status stat[2]; 
  MPI_Datatype MPI_ROW_BFR;
  MPI_Type_vector(0.5*(par_rows - 2), 1, 1, MPI_DOUBLE, &MPI_ROW_BFR);
  MPI_Type_commit(&MPI_ROW_BFR);

  if(crds[1] == Q/2 -1 && crds[0] < P/2) {

    int first_target = refine_ranks[crds[0]*2];
    int second_target = refine_ranks[(crds[0]*2)+1];
    // send first half
    MPI_Send(&mat[par_rows - 2][1], 1, MPI_ROW_BFR, first_target, 0, MPI_COMM_WORLD);
    
    // send second half
    MPI_Send(&mat[par_rows - 2][par_rows/2], 1, MPI_ROW_BFR, second_target, 0, MPI_COMM_WORLD);
  }

  if(rcrds[1] == 0) {
    MPI_Recv(&buffer[0], 1, MPI_ROW_BFR, coarse_ranks[ P*((Q/2) -1) + rcrds[0]/2], 0, MPI_COMM_WORLD, &stat[1]);
  }

}

//sends right interface points from refined grid to coarse grid
void sendRightInterfacePoints(double** mat, double* buffer, int par_ref_rows, int par_ref_cols)
{
  MPI_Status stat[2]; 
  MPI_Datatype MPI_CLM_BFR;
  MPI_Type_vector(par_ref_rows - 2, 1, par_ref_cols, MPI_DOUBLE, &MPI_CLM_BFR);
  MPI_Type_commit(&MPI_CLM_BFR);
/*
  if(rcrds[0] == P-1) 
  {
    if (    )
      MPI_Send(&mat[1][par_ref_cols - 2], 1, MPI_CLM_BFR, 3, 0, MPI_COMM_WORLD);
    else
      MPI_Send(&mat[1][par_ref_cols - 2], 1, MPI_CLM_BFR, 3, 0, MPI_COMM_WORLD);

  }
*/
  if(crds[0] == P/2 && crds[1] >= Q/2) 
  {
  MPI_Recv(&buffer[0], par_ref_rows-2, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD, &stat[1]);
  }

}

//sends top interface points from refined grid to coarse grid
void sendTopInterfacePoints(double** mat, double* buffer, int par_ref_rows, int par_ref_cols)
 {
  MPI_Status stat[2]; 
  MPI_Datatype MPI_ROW_BFR;
  MPI_Type_vector(par_ref_cols - 2, 1, 1, MPI_DOUBLE, &MPI_ROW_BFR);
  MPI_Type_commit(&MPI_ROW_BFR);

/*
  if(rcrds[1] == 0)
  {
    if (    )
      MPI_Send(&mat[1][1], 1, MPI_ROW_BFR, coarse_ranks[ P*((Q/2) -1) + rcrds[0]/2], 0, MPI_COMM_WORLD);
    else
      MPI_Send(&mat[1][0], 1, MPI_ROW_BFR, coarse_ranks[ P*((Q/2) -1) + rcrds[0]/2], 0, MPI_COMM_WORLD);
  }
*/
  if(crds[1] == Q/2 -1 && crds[0] < P/2) 
  {
    int first_target = refine_ranks[crds[0]*2];
    int second_target = refine_ranks[(crds[0]*2)+1];
    
    MPI_Recv(&buffer[0], 1, MPI_ROW_BFR, first_target, 0, MPI_COMM_WORLD, &stat[1]);
    MPI_Recv(&buffer[0], 1, MPI_ROW_BFR, second_target, 0, MPI_COMM_WORLD, &stat[2]);
  }
}

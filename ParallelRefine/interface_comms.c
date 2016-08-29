#include "interface_comms.h"
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
  if(crds[0] == Q/2 -1 && crds[1] < P/2) {
    for (i = 1; i < par_cols - 1; i++)
      {
	mat[par_rows - 1][i] = toparray[j];
	j++;
      }
  }

  j=0;
  //right
  if(crds[1] == P/2 && crds[0] >= Q/2) {
    for(i = 1; i < par_rows - 1; i++)
      {
	mat[i][0] = rightarray[j];
	j++;
      }
  }
}

//injects buffers into refined grid
void injectFinePoints(double* topbuffer, double* rightbuffer, double** mat, int par_ref_rows, int par_ref_cols)
{
  int i, j = 0, n = 0;

  int start_top=1;
  int start_right=1;	

  if(par_ref_cols%2 == 1 && rcrds[1]%2 == 0)
    start_top = 0; // if we start at an odd boundary

  if(par_ref_rows%2 == 1 && rcrds[0]%2 == 1)
    start_right = 0; // if we start at an odd boundary

  //  if(par_ref_cols%2 == 0 && rcrds[1]%2 == 1)    start_top = 2;
  if(par_ref_cols%2 == 0) 
    {
      start_top = 0;
      /* if(rcrds[1]%2 == 0) */
      /* 	{ */
      /* 	  start_top = 0; */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  start_top = 2; */
      /* 	} */
    }

  //top
  if(rcrds[0] == 0) {
    for (i = start_top; i < par_ref_cols; i+=2)
      {
	mat[0][i] = topbuffer[j];
	j++;
      }

    //average in between points
    for(i=start_top+1; i<par_ref_cols-1; i+=2)
      {
    	mat[1][i] = 0.5*(mat[1][i-1] + mat[1][i+1]);
      }
  }

  if(rcrds[1] == P-1) {

    for(i = start_right; i < par_ref_rows; i+=2)
      {
	mat[i][par_ref_cols - 1] = rightbuffer[n];
	n++;
      }
    
    // average in between points
    for(i=start_right+1; i<par_ref_rows-1; i+=2)
      {
	mat[i][par_ref_cols-2] = 0.5*(mat[i-1][par_ref_cols-2] + mat[i+1][par_ref_cols-2]);
      }
  }
}

//sends points from coarse grid to right of refined grid
void sendRightPoints(double** mat, double* buffer, int par_rows, int par_cols)
{
  int i;
  MPI_Status stat; 
  MPI_Datatype MPI_CLM_BFR;
  double* first_bfr = calloc(par_rows/2, sizeof(double));
  double* second_bfr = calloc(par_rows/2, sizeof(double));
  MPI_Type_vector(par_rows/2, 1, 2, MPI_DOUBLE, &MPI_CLM_BFR);
  MPI_Type_commit(&MPI_CLM_BFR);

  for(i = 0; i < par_rows/2; i++) {
    first_bfr[i] = mat[i+1][1];
    second_bfr[i] = mat[par_rows/2 + i][1];
  }

  if(crds[1] == P/2 && crds[0] >= Q/2) {
  
    int first_target  = refine_ranks[  P-1 + 2*P*(crds[0]-Q/2) ];
    int second_target = refine_ranks[2*P-1 + 2*P*(crds[0]-Q/2) ];
    
    MPI_Send(first_bfr,  par_rows/2, MPI_DOUBLE, first_target, 0, MPI_COMM_WORLD);
    MPI_Send(second_bfr, par_rows/2, MPI_DOUBLE, second_target, 0, MPI_COMM_WORLD);
  }

  if(rcrds[1] == P-1) {
   
    MPI_Recv(buffer, par_rows/2, MPI_DOUBLE, coarse_ranks[(Q/2)*P + (P/2)*(rcrds[0]/2) ], 0, MPI_COMM_WORLD, &stat);
  }

}


//sends points from coarse grid to top of refined grid
void sendTopPoints(double** mat, double* buffer, int par_rows, int par_cols)
{
  MPI_Status stat; 

  if(crds[0] == Q/2 -1 && crds[1] < P/2) {

    int first_target  = refine_ranks[ crds[1]*2];
    int second_target = refine_ranks[(crds[1]*2)+1];
    
    // send first half
    MPI_Send(&mat[par_rows - 2][0], par_cols/2, MPI_DOUBLE, first_target, 0, MPI_COMM_WORLD);
    
    // send second half 
    // CHECK THIS, CHANGED IN MY SLEEP -- before par_cols/2
    MPI_Send(&mat[par_rows - 2][/**/(par_cols-1)/2/**/], par_cols/2, MPI_DOUBLE, second_target, 0, MPI_COMM_WORLD);
  }

  if(rcrds[0] == 0) {
    MPI_Recv(buffer, par_cols/2, MPI_DOUBLE, coarse_ranks[ P*((Q/2) -1) + rcrds[1]/2], 0, MPI_COMM_WORLD, &stat);
  }

}


//sends right interface points from refined grid to coarse grid
void sendRightInterfacePoints(double** mat, double* buffer, int par_ref_rows, int par_ref_cols)
{
  int i,j;
  MPI_Status stat[2]; 
  double* bfr = calloc((par_ref_rows-1)/2, sizeof(double));

  if(rcrds[1] == P-1) 
  {
    int start;
    if (rcrds[0]%2 == 0 || par_ref_rows%2 == 0)
      {
	start = 1;
      }
    else
      {
	start = 2;
      }

    j=0;
    for(i=start; i<par_ref_rows-1; i+=2) 
      {
	bfr[j] = mat[i][par_ref_cols-2];
	j++;
      }
  
        MPI_Send(bfr, (par_ref_rows-1)/2, MPI_DOUBLE, coarse_ranks[(Q/2)*P + (P/2)*(rcrds[0]/2) ], 0, MPI_COMM_WORLD);
  }
 
  if(crds[1] == P/2 && crds[0] >= Q/2)
  {
    int first_target  = refine_ranks[  P-1 + 2*P*(crds[0]-Q/2) ];
    int second_target = refine_ranks[2*P-1 + 2*P*(crds[0]-Q/2) ];

    MPI_Recv(buffer,                      (par_ref_rows-1)/2, MPI_DOUBLE, first_target,  0, MPI_COMM_WORLD, &stat[0]);
    MPI_Recv(buffer + (par_ref_rows-1)/2, (par_ref_rows-1)/2, MPI_DOUBLE, second_target, 0, MPI_COMM_WORLD, &stat[1]);
  }

}

//sends top interface points from refined grid to coarse grid
void sendTopInterfacePoints(double** mat, double* buffer, int par_ref_rows, int par_ref_cols)
 {
   int i,j;
   MPI_Status stat[2];
   double* bfr = calloc((par_ref_cols-1)/2, sizeof(double)); 

   if(rcrds[0] == 0) 
     {
       int start;
       if (rcrds[1]%2 == 0 || par_ref_rows%2 == 0)
	 {
	   start = 2;
	 }
       else
	 {
	   start = 1;
	 }
       
       i=0;
       for(j=start; j<par_ref_cols-1; j+=2) 
	 {
	   bfr[i] = mat[1][j];
	   i++;
	 }
     }


  if(rcrds[0] == 0)
    {                                                           // P*( (Q/2) -1 ) + rcrds[1]/2
    MPI_Send(bfr, (par_ref_cols-1)/2, MPI_DOUBLE, coarse_ranks[P*( (Q/2) -1 ) + rcrds[1]/2], 0, MPI_COMM_WORLD);
  }

  if(crds[0] == Q/2 -1 && crds[1] < P/2) 
  {
    int first_target =  refine_ranks[crds[1]*2]; // changed from crds[0] to crds[1]
    int second_target = refine_ranks[(crds[1]*2)+1]; // same
    
    MPI_Recv(buffer,                        (par_ref_cols-1)/2, MPI_DOUBLE, first_target,  0, MPI_COMM_WORLD, &stat[0]);
    MPI_Recv(buffer + (par_ref_cols - 2)/2, (par_ref_cols-1)/2, MPI_DOUBLE, second_target, 0, MPI_COMM_WORLD, &stat[1]);
  }
}

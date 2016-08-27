
#include "compute_timestep.h"
#include <mpi.h>
#include <math.h>

extern int P, Q, rank, size;
extern int crds[2], rcrds[2];
extern MPI_Comm MPI_COMM_CART, MPI_COMM_CART2;
extern MPI_Datatype MPI_CLM, MPI_CLM_FINE;

/****    Interface Points *****/
void computeCornerTimestep(double* buffer, double** mat1, double** mat2, double wxy, int par_ref_cols)
{
  int i = 1, j = par_ref_cols - 2;

  mat2[i][j] = /*prev pointmat1_corner*/mat1[i][j] + wxy * (16/15) * (/*prev point*/-4* mat1[i][j]/*mat1_corner*/+
		.5 * /*bound lwr*/mat1[i+2][j] + /*coarse*/mat1[i][j+1] + /*bound*/ mat1[i-1][j] + .5 * /*coarse*/mat1[i][j-2]
	        +/*mesh*/ mat1[i+1][j-1]);

}



void computeInterfaceRightTimestep(double* buffer, double** mat1, double** mat2, int par_ref_rows, int ncols, double wxy, double dx, int refcols)
{
  int j = par_ref_rows - 2, i;

  int u,d;
  double upval=0.0, downval=0.0;

  MPI_Status stat;
  MPI_Cart_shift(MPI_COMM_CART2, 1, 1, &u, &d);

  int start = 1;

  if(rcrds[1] == Q/2)
    start = 3;
  else if(par_ref_rows%2 == 1 && rcrds[1]%2 == 1)
    start = 2;

  // send/recv interface values

  // if cols even
  // recv one value from down
  // send one value to up

  // if cols odd
      // if coord even
      // send/recv to both up and down
      
      // if coord odd
      // do nothing
  

  if(par_ref_rows%2 == 0)
    {
      upval = mat1[0][j];
      MPI_Sendrecv(&mat1[2][j], 1, MPI_DOUBLE, u, 0, &downval, 1, MPI_DOUBLE, d, 0, MPI_COMM_CART2, &stat);
    }
  else
    {
      if(rcrds[0]%2 == 0)
	{
	  MPI_Sendrecv(&mat1[2][j],              1, MPI_DOUBLE, u, 0, &downval, 1, MPI_DOUBLE, d, 0, MPI_COMM_CART2, &stat);
	  MPI_Sendrecv(&mat1[par_ref_rows-3][j], 1, MPI_DOUBLE, d, 0, &upval,   1, MPI_DOUBLE, u, 0, MPI_COMM_CART2, &stat);
	}
      else
	{
	  upval =   mat1[0][j];
	  downval = mat1[par_ref_rows-1][j];
	}
    }

  // calculate first interface point

      mat2[start][j] =/*previous interface point*/ mat1[start][j] + /*previous interface point*/ mat1[start][j] * wxy * (dx/2) + wxy *
	(4/3) * (/*previous interface point*/- 4*mat1[start][j] + /*coarse*/mat1[start][j+1]  + /*coarse*/ .5 * mat1[start+2][j]
		 +/*coarse*/ .5 * upval +/*refine*/ .5 * mat1[start-1][j-1] +/*refine*/ mat1[start][j-1] +/*refine*/ .5
		 * mat1[start+1][j-1]);


  for (i = start+2; i < par_ref_rows - 2; i+=2)
    {
      mat2[i][j] =/*previous interface point*/ mat1[i][j] + /*previous interface point*/ mat1[i][j] * wxy * (dx/2) + wxy *
	(4/3) * (/*previous interface point*/- 4*mat1[i][j] + /*coarse*/mat1[i][j+1]  + /*coarse*/ .5 * mat1[i+2][j]
		 +/*coarse*/ .5 * mat1[i-2][j] +/*refine*/ .5 * mat1[i-1][j-1] +/*refine*/ mat1[i][j-1] +/*refine*/ .5
		 * mat1[i+1][j-1]);
    }

  // dont update extremal interface value - its a constant boundary
  if( rcrds[1] != Q-1 ) {
      mat2[par_ref_rows-2][j] =/*previous interface point*/ mat1[par_ref_rows-2][j] + /*previous interface point*/ mat1[par_ref_rows-2][j]
      * wxy * (dx/2) + wxy * (4/3) * (/*previous interface point*/- 4*mat1[par_ref_rows-2][j] + /*coarse*/mat1[par_ref_rows-2][j+1]  
      + /*coarse*/ .5 * downval +/*coarse*/ .5 * mat1[par_ref_rows-4][j] +/*refine*/ .5 * mat1[par_ref_rows-3][j-1] 
      +/*refine*/ mat1[par_ref_rows-2][j-1] +/*refine*/ .5 * mat1[par_ref_rows-1][j-1]);
  }

}


void computeInterfaceTopTimestep(double* buffer, double** mat1, double** mat2, int par_rows, int par_ref_cols, double wxy, double dx)
{
  int i = 1;
  int j;
  int l,r;
  double leftval, rightval;

  MPI_Status stat;
  MPI_Cart_shift(MPI_COMM_CART2, 1, 1, &l, &r);

  int start = 2;

  if(rcrds[0] == 0)
    start = 4;
  else if(par_ref_cols%2 == 1 && rcrds[0]%2 == 1)
    start = 1;

  // send/recv interface values

  // if cols even
  // recv one value from the right
  // send one value to the left

  // if cols odd
      // if coord odd
      // send/recv to both left and right
      
      // if coord even
      // do nothing
  

  if(par_ref_cols%2 == 0)
    {
      leftval = mat1[i][0];
      MPI_Sendrecv(&mat1[1][2], 1, MPI_DOUBLE, l, 0, &rightval, 1, MPI_DOUBLE, r, 0, MPI_COMM_CART2, &stat);
    }
  else
    {
      if(rcrds[0]%2 == 1)
	{
	  MPI_Sendrecv(&mat1[1][2],              1, MPI_DOUBLE, l, 0, &rightval, 1, MPI_DOUBLE, r, 0, MPI_COMM_CART2, &stat);
	  MPI_Sendrecv(&mat1[1][par_ref_cols-3], 1, MPI_DOUBLE, r, 0, &leftval,  1, MPI_DOUBLE, l, 0, MPI_COMM_CART2, &stat);
	}
      else
	{
	  leftval = mat1[i][0];
	  rightval = mat1[i][par_ref_cols-1];
	}
    }

  // calculate first interface point

  mat2[i][start] = /*previous interface point*/ mat1[i][start] + /*prev interface point*/ mat1[i][start] * wxy * (dx/2) + wxy * (4/3) *
    (-/*prev interface point*/ 4*mat1[i][start] +/*coarse*/ mat1[i-1][start]  +/*coarse*/ .5 * mat1[i][start+2]
     +/*coarse*/ .5 * leftval +/*refine*/ .5 * mat1[i+1][start-1] +/*refine*/ mat1[i+1][start] +/*refine*/ .5 * mat1[i+1][start+1]);


  // calculate intermediary interface points
  for (j = start+2; j < par_ref_cols - 2; j+=2)
    {
        mat2[i][j] = /*previous interface point*/ mat1[i][j] + /*prev interface point*/ mat1[i][j] * wxy * (dx/2) + wxy * (4/3) *
	(-/*prev interface point*/ 4*mat1[i][j] +/*coarse*/ mat1[i-1][j]  +/*coarse*/ .5 * mat1[i][j+2]
	 +/*coarse*/ .5 * mat1[i][j-2] +/*refine*/ .5 * mat1[i+1][j-1] +/*refine*/ mat1[i+1][j] +/*refine*/ .5 * mat1[i+1][j+1]);
    }

  // if we're not at the corner

  if(rcrds[0] != P-1) {
    // calcuate last interface point
    mat2[i][par_ref_cols-2] = /*previous interface point*/ mat1[i][par_ref_cols-2] + /*prev interface point*/ mat1[i][par_ref_cols-2] 
    * wxy * (dx/2) + wxy * (4/3) * (-/*prev interface point*/ 4*mat1[i][par_ref_cols-2] +/*coarse*/ mat1[i-1][par_ref_cols-2]  
    +/*coarse*/ .5 * rightval +/*coarse*/ .5 * mat1[i][par_ref_cols-4] +/*refine*/ .5 * mat1[i+1][par_ref_cols-3] 
    +/*refine*/ mat1[i+1][par_ref_cols-2] +/*refine*/ .5 * mat1[i+1][par_ref_cols-1]);
  }
}

void computeFineTimestep(double** mat1_refine, double** mat2_refine, int par_ref_rows, int par_ref_cols, double wxy, double sxy)
{
  int i,j;
/*
  MPI_Status stat[8]; 
  MPI_Request req[8];
  int u, d, l, r;
  int startrow, startcol, endrow, endcol;

  MPI_Cart_shift(MPI_COMM_CART2, 0, 1, &u, &d);
  MPI_Cart_shift(MPI_COMM_CART2, 1, 1, &l, &r);

  startcol = 1;
  endcol = par_ref_cols - 1;
  startrow = 1;
  endrow = par_ref_rows - 1;

  if (rcrds[1] == 0)
    startcol = 2;
  else if (rcrds[1] == S - 1)
    endcol = par_ref_cols -2;

  if (rcrds[0] == 0)
    startrow = 2;
  else if (rcrds[0] == R - 1)
    endrow = par_ref_rows -2;

  //solve edges
  if (rcrds[0] != 0)
    {
      //not in top blocks
      i = 1;
      for (j = startcol; j < endcol; j++)
	{
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
	     -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	    }
	}
      
  if (rcrds[0] != R-1)
    {
      //not in bottom blocks
      i = par_ref_rows - 2;
      for (j = startcol; j < endcol; j++)
	{
	  mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
    -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	}
    }

  if (rcrds[1] != 0)
    {
      j = 1;
      for (i = startrow; i < endrow; i++)
	{
	  mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
     -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	}
    }


  if (rcrds[1] != S-1)
    {
      j = par_ref_cols - 2;
      for (i = startrow; i < endrow; i++)
	{
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
	   -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	}
    }

  //halo exchange
  // Send rows up and down
  MPI_Irecv(mat2_refine[0],              par_ref_cols, MPI_DOUBLE, u, 0, MPI_COMM_CART2, &req2[0]);
  MPI_Irecv(mat2_refine[par_ref_rows-1], par_ref_cols, MPI_DOUBLE, d, 1, MPI_COMM_CART2, &req2[1]);
  MPI_Isend(mat2_refine[par_ref_rows-2], par_ref_cols, MPI_DOUBLE, d, 0, MPI_COMM_CART2, &req2[2]);
  MPI_Isend(mat2_refine[1],              par_ref_cols, MPI_DOUBLE, u, 1, MPI_COMM_CART2, &req2[3]);
  // Send columns left and right
  MPI_Irecv(&mat2_refine[1][par_ref_cols-1], 1, MPI_CLM, r, 2, MPI_COMM_CART2, &req2[4]);
  MPI_Irecv(&mat2_refine[1][0],              1, MPI_CLM, l, 3, MPI_COMM_CART2, &req2[5]);
  MPI_Isend(&mat2_refine[1][1],              1, MPI_CLM, l, 2, MPI_COMM_CART2, &req2[6]);
  MPI_Isend(&mat2_refine[1][par_ref_cols-2], 1, MPI_CLM, r, 3, MPI_COMM_CART2, &req2[7]);

  for (i = startrow; i < endrow; i++)
    {
      for (j = startcol; j < endcol; j++)
	{
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
	     -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	}
    }

  MPI_Waitall(8, req2, stat2);
*/

  for (i = 2; i < par_ref_rows; i++)
      {
      for (j = 2; j < par_ref_cols; j++)
	{
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
	      -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	}
    }
}


void computeTimestep(double* bfr_top_refine_pts, double* bfr_right_refine_pts, double** mat1, double** mat2, int par_rows, int par_cols, double wx, double wy)
{
  int i, j;
  MPI_Status stat[8]; 
  MPI_Request req[8];
  int u, d, l, r;
  int startrow, startcol, endrow, endcol;


  MPI_Cart_shift(MPI_COMM_CART, 0, 1, &u, &d);
  MPI_Cart_shift(MPI_COMM_CART, 1, 1, &l, &r);

  
  startcol = 1;
  endcol = par_cols - 1;
  startrow = 1;
  endrow = par_rows - 1;

  //sets where on grid to start calculations
  if (crds[1] == 0)
    startcol = 2;
  else if (crds[1] == P - 1)
    endcol = par_cols -2;

  if (crds[0] == 0)
    startrow = 2;
  else if (crds[0] == Q - 1)
    endrow = par_rows -2;

 
  //solve edges
  if (crds[0] != 0)
    {
      //not in top blocks
      i = 1;
      for (j = startcol; j < endcol; j++)
	{
	  mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
	}
    }

  if (crds[0] != Q-1)
    {
      //not in bottom blocks
      i = par_rows - 2;
      for (j = startcol; j < endcol; j++)
	{
	  mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
	}
    }

  if (crds[1] != 0)
    {
      //not in left blocks
      j = 1;
      for (i = startrow; i < endrow; i++)
	{
	  mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
	}
    }

  if (crds[1] != P-1)
    {
      //not in right blocks
      j = par_cols - 2;
      for (i = startrow; i < endrow; i++)
	{
	  mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
	}
    }
  
  //halo exchange
  // Send rows up and down
  MPI_Irecv(mat2[0],      	par_cols, MPI_DOUBLE, u, 0, MPI_COMM_CART, &req[0]);
  MPI_Irecv(mat2[par_rows-1], 	par_cols, MPI_DOUBLE, d, 1, MPI_COMM_CART, &req[1]);
  MPI_Isend(mat2[par_rows-2], 	par_cols, MPI_DOUBLE, d, 0, MPI_COMM_CART, &req[2]);
  MPI_Isend(mat2[1],      	par_cols, MPI_DOUBLE, u, 1, MPI_COMM_CART, &req[3]);
  // Send columns left and right
  MPI_Irecv(&mat2[1][par_cols-1],	 1, MPI_CLM, r, 2, MPI_COMM_CART, &req[4]);
  MPI_Irecv(&mat2[1][0],      		 1, MPI_CLM, l, 3, MPI_COMM_CART, &req[5]);
  MPI_Isend(&mat2[1][1],      		 1, MPI_CLM, l, 2, MPI_COMM_CART, &req[6]);
  MPI_Isend(&mat2[1][par_cols-2],	 1, MPI_CLM, r, 3, MPI_COMM_CART, &req[7]);

  //interior calculations
  if(rank != 2){
    for (i = startrow; i < endrow; i++)
      {
	for (j = startcol; j < endcol; j++)
	  {
	    mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
	  }}
  }
  //ranks in bottom row must calculate from halo row
  //if (rank == 0) 
  if(crds[1] == Q/2 -1 && crds[0] < P/2)
    {
      for(j = 2; j < par_cols - 1; j++)
      {
	mat2[par_rows-2][j] = mat1[par_rows-2][j] + wx * (mat1[par_rows-2+1][j] - 2*mat1[par_rows-2][j] + mat1[par_rows-2-1][j]) + wy * (mat1[par_rows -2][j+1] - 2 * mat1[par_rows-2][j] + mat1[par_rows-2][j-1]);
      }
    }

  //ranks in right row must calculate using halo row
  //if (rank == 3){
  if(crds[0] == P/2 && crds[1] >= Q/2)
  {
    for( i = 1; i < par_rows - 2; i++)
      mat2[i][1] = mat1[i][1] + wx * (mat1[i+1][1] - 2*mat1[i][1] + mat1[i-1][1]) + wy * (mat1[i][1+1] - 2 * mat1[i][1] + mat1[i][1-1]);
  }

  MPI_Waitall(8,req, stat);
}
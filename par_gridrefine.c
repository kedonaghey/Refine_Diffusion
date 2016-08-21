#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <mpi.h>

int rank, size, cRank, cRank2;
MPI_Comm MPI_COMM_CART;
MPI_Comm MPI_COMM_CART2;
int P, Q, R, S;
int crds[2];
int rcrds[2];
MPI_Datatype MPI_CLM;
MPI_Datatype MPI_CLM_FINE;

void initFineGrid(double** mat, int nrows, int ncols)
{
  int i;
  int bottom = 50;
  int left = 20;
  // left rows initial conditions
  if(rcrds[1] == 0)
    {
      for ( i = 1; i < nrows - 1; i+=2)
	{
	  mat[i][1] = left;
	}
      for (i = 2; i < nrows - 1; i+=2)
	{
	  mat[i][1] = left/2;
	}
    }

  // bottom rows initial conditions
  if(rcrds[0] == S-1)
    {
      for ( i = 1; i < ncols - 1; i+=2)
	{
	  mat[nrows-2][i] = bottom;//100
	}
      for (i = 2; i < ncols - 1; i+=2)
	{
	  mat[nrows-2][i] = bottom/2;
	}
    } 
}

void initGrid(double** mat, int nrows, int ncols, int crse_rows_cutoff, int crse_cols_cutoff)
{
  int i;
  //left and right rows initial conditions
  if(crds[1] == 0)
    {
      //  for ( i = 1; i < nrows - crse_rows_cutoff +1; i++)
      for ( i = 1; i < nrows - 1; i++)
	{
	  mat[i][1] = 20;//0
	}
    }
  else if(crds[1] == P - 1)
    {
      for (i = 1; i <nrows - 1; i++)
	{ 
	  mat[i][ncols - 2] = 30;//100
	}
    }

  //bottom and top rows initial conditions
  if(crds[0] == 0)
    {
      for ( i = 1; i < ncols - 1; i++)
	{
	  mat[1][i] = 50;//0
	}
    }
  else if(crds[0] == Q-1)
    {
      for ( i = 1; i < ncols - 1; i++)
	{
	  mat[nrows-2][i] = 40;//100
	}
    }

}

void printGrid(double** mat, int nrows, int ncols)
{
  int i, j, k;

  for(k=0; k<size; k++) {
    if(rank==k) {
      printf("\nRank %d:\n", rank);
      /* printf("%f\n ", mat[1][1]); */
      /* printf("%f\n ", mat[5][cols-2]); */
      for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
          printf("%f ", mat[i][j]);
        }
        printf("\n");
      }
    }
    usleep(200);
  }


}

void injectCoarsePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, int crse_rows_cutoff, int crse_cols_cutoff, double corner)
{
  int i, j = 0;
  //int test[] = {1,2,3};
  //top
  for (i = 1; i < crse_cols_cutoff; i++)
    mat[nrows-crse_rows_cutoff][i] = toparray[i-1];//toparray[i-1];//0
  //right
  for(i = nrows-2; i > nrows - crse_rows_cutoff; i--)
    {
      //mat[i][crse_cols_cutoff - 1] = test[i-crse_rows_cutoff-1];//rightarray[i-crs-1];//100
      mat[i][crse_cols_cutoff - 1] = rightarray[j];
      j++;
    }
  //corner
  mat[nrows-crse_rows_cutoff][crse_cols_cutoff-1] = corner;
}

void injectFinePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, double corner)
{
  int i, j = 0, n = 0;
  //top
  //int test[] = {1,2,3};
  for (i = 2; i < ncols - 1; i+=2)
    {
      mat[0][i] = toparray[j];//toparray[i-1];//0
      j++;
    }
  //right
  for(i = nrows - 3; i > 0; i-=2)
    {
      //mat[i][ncols - 1] = rightarray[n];//rightarray[i-1];//100
      mat[i][ncols - 1] = rightarray[n];
      n++;
    }
  //corner
  //double check
  mat[0][ncols-1] = corner;
}

/****    Interface Points *****/
//for bound points use regular matrix
//finished
double computeCornerTimestep(double** mat1_refine, double** mat1, int i, int j, int m, int n, double wxy)
{
  //double mat1_corner;
  double mat2_corner;

  mat2_corner = /*prev pointmat1_corner*/mat1[i][j] + wxy * (16/15) * (/*prev point*/-4* mat1[i][j]/*mat1_corner*/+
								       .5 * /*bound lwr*/mat1[i+1][j] + /*coarse*/mat1[i][j+1] + /*bound*/ mat1[i-1][j] + .5 * /*coarse*/mat1[i][j-1]
								       +/*mesh*/ mat1_refine[m+1][n-1]);

  return mat2_corner;
}

//need to add in different

void computeInterfaceRightTimestep(double* rightarray, double** mat1_refine, double** mat1, int nrows, int ncols, int crse_cols_cutoff, int m, double wxy, double dx, int crse_rows_cutoff, int refcols)
{
  int j = crse_cols_cutoff-1, n = refcols-1, t = 0, i;
  m = 2;

  for (i = nrows-2; i > nrows-crse_rows_cutoff; i--)
    {
      rightarray[t] =/*previous interface point*/ mat1[i][j] + /*previous interface point*/ mat1[i][j] * wxy * (dx/2) + wxy *
	(4/3) * (/*previous interface point*/- 4*mat1[i][j] + /*coarse*/mat1[i][j+1]  + /*coarse*/ .5 * mat1[i+1][j]
		 +/*coarse*/ .5 * mat1[i-1][j] +/*refine*/ .5 * mat1_refine[m-1][n-1] +/*refine*/ mat1_refine[m][n-1] +/*refine*/ .5
		 * mat1_refine[m+1][n-1]);//changed one n to - different from the book
      m+=2;
      //iterates through right array
      t++;
    }
}


void computeInterfaceTopTimestep(double* toparray, double** mat1_refine, double** mat1, int nrows, int ncols, int crse_rows_cutoff, int n, double wxy, double dx, int crse_cols_cutoff)
{
  int i = nrows - crse_rows_cutoff;
  int m = 0, t = 0, j;
  n = 2;

  for (j = 1; j < crse_cols_cutoff; j++)
    {
      toparray[t] = /*previous interface point*/mat1[i][j] + /*prev interface point*/mat1[i][j] * wxy * (dx/2) + wxy * (4/3) * \
	(-/*prev interface point*/ 4*mat1_refine[i][j] +/*coarse*/ mat1[i-1][j]  +/*coarse*/ .5 * mat1_refine[i][j+1]\
	 +/*coarse*/ .5 * mat1_refine[i][j-1] +/*refine*/ .5 * mat1_refine[m+1][n-1] +/*refine*/ mat1_refine[m+1][n] +/*refine*/ .5 * mat1_refine[m+1][n+1]);
      n+=2;
      //iterates through top array
      t++;
    }

}


void sendRightRefinePoints(double** mat, double* buffer, int par_ref_rows, int par_ref_cols)
{
  MPI_Status stat[2]; 
  MPI_Request req[2];
  MPI_Datatype MPI_CLM_BFR;
  MPI_Type_vector(par_ref_rows-2, 1, par_ref_cols, MPI_DOUBLE, &MPI_CLM_BFR);
  MPI_Type_commit(&MPI_CLM_BFR);

  if(rank==4)
  MPI_Send(&mat[1][par_ref_cols - 3], 1, MPI_CLM_BFR, 3, 0, MPI_COMM_WORLD);

  if(rank==3)
  MPI_Recv(&buffer[0], par_ref_rows-2, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD, &stat[1]);

  printf("%d in refine\n", rank);

  //MPI_Waitall(2, req, stat);
}

void sendTopRefinePoints(double** mat, double* buffer, int par_ref_rows, int par_ref_cols)
{
  MPI_Status stat[2]; 
  MPI_Request req[2];
/*MPI_Datatype MPI_CLM_BFR;
  MPI_Type_vector(par_ref_rows-2, 1, par_ref_cols, MPI_DOUBLE, &MPI_CLM_BFR);
  MPI_Type_commit(&MPI_CLM_BFR);
*/
  if(rank==4)
  MPI_Send(mat[1], par_ref_cols , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

  if(rank==0)
  MPI_Recv(&buffer[0], par_ref_cols, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD, &stat[1]);

  printf("%d in refine\n", rank);

  //MPI_Waitall(2, req, stat);
}

void computeFineTimestep(double** mat1_refine, double** mat2_refine, int par_ref_rows, int par_ref_cols, double* converge, double wxy, double sxy)
{

  int i,j;
  MPI_Status stat2[8]; 
  MPI_Request req2[8];
  int u, d, l, r;
  int startrow, startcol, endrow, endcol;
/*
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
	  if (j%2 == 1)
	    {
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i-1][j+1]
	     + mat1_refine[i-1][j-1]) -2 * mat1_refine[i][j] + mat1_refine[i+1][j]) + wxy * (mat1_refine[i][j-1] -2 * mat1_refine[i][j]																	 + mat1_refine[i][j+1]);
	    }
	  else
	    {
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
	     -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	    }
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
	  if (i%2 == 1)
	    {
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i+1][j+1]
	  + mat1_refine[i-1][j+1]) -2 * mat1_refine[i][j] + mat1_refine[i][j-1]) + wxy * (mat1_refine[i+1][j] -2 * mat1_refine[i][j]
	  + mat1_refine[i-1][j]);
	    }
	  else
	    {
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
	   -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	    }
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
	  if (i == 2 && j%2 == 0)
	    {
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i-1][j+1]
	      + mat1_refine[i-1][j-1]) -2 * mat1_refine[i][j] + mat1_refine[i+1][j]) + wxy * (mat1_refine[i][j-1] -2 * mat1_refine[i][j]
	      + mat1_refine[i][j+1]);
	    }
	  if (j == par_ref_cols-3 && i%2 == 0)
	    {
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i+1][j+1]
	      + mat1_refine[i-1][j+1]) -2 * mat1_refine[i][j] + mat1_refine[i][j-1]) + wxy * (mat1_refine[i+1][j] -2 * mat1_refine[i][j]\
	      + mat1_refine[i-1][j]);
	    }
	  else
	    {
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
	     -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	    }
	}
    }

  MPI_Waitall(8, req2, stat2);
*/

  for (i = 2; i < par_ref_rows; i++)
      {
      for (j = 2; j < par_ref_cols; j++)
	{

/*	  if (i == 2 && j%2 == 0)
	    {
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i-1][j+1]
	      + mat1_refine[i-1][j-1]) -2 * mat1_refine[i][j] + mat1_refine[i+1][j]) + wxy * (mat1_refine[i][j-1] -2 * mat1_refine[i][j]
	      + mat1_refine[i][j+1]);
	    }
	  if (j == par_ref_cols-3 && i%2 == 0)
	    {
	      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i+1][j+1]
	      + mat1_refine[i-1][j+1]) -2 * mat1_refine[i][j] + mat1_refine[i][j-1]) + wxy * (mat1_refine[i+1][j] -2 * mat1_refine[i][j]\
	      + mat1_refine[i-1][j]);
	    }
	  else
*/	    {
	      mat2_refine[i][j] = 1;//mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]
	      //-2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
	    }
	}
    }
}


void computeTimestep(double** mat1, double** mat2, int par_rows, int par_cols, double* converge, double wx, double wy, int crse_rows_cutoff, int crse_cols_cutoff)
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
      j = 1;
      for (i = startrow; i < endrow; i++)
	{
	  mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
	}
    }

  if (crds[1] != P-1)
    {
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


  for (i = startrow; i < endrow; i++)
  {
  for (j = startcol; j < endcol; j++)
  {
  mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
  }}

  MPI_Waitall(8,req, stat);
  
}

void update(double*** mat1_ptr, double*** mat2_ptr)
{
  double **tmp = *mat1_ptr;
  *mat1_ptr = *mat2_ptr;
  *mat2_ptr = tmp;
}

int main(int argc, char *argv[])
{

  int i, c;
  double **mat1, **mat2;
  //  double **mat1_refine, **mat2_refine;
  //size of grid
  int nrows = 10, ncols = 10;
  //2x2 for mesh in a 4x4 array
  double dx = 0, dy= 0, dt, time;
  double converge = 0;
  int iter, max_iter = 190;
  //clock_t start, end;
  //int prnt = 1;
  char output_file[80] = "k.txt";
  FILE *fp;

  //MPI Init
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  P = 2;
  Q = 2;

  R = 1;
  S = 1;

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

  //parallel rows & cols coarse grid
  int par_rows = 2 + nrows/Q;
  int par_cols = 2 + ncols/P;

  int szofmesh = .5*nrows;

  //size of fine grid
  int refrows = szofmesh * 2 + 1;
  int refcols = szofmesh * 2 + 1;

  //parallel refined rows & cols
  int par_ref_rows = 2 + refrows/R;
  int par_ref_cols = 2 + refcols/S;


  //  MPI_Type_vector(par_crows - 2, 1, par_cols * 2, MPI_DOUBLE, &MPI_CLM);
  //  MPI_Type_commit(&MPI_ROW);

  MPI_Group coarse_grid_group;
  MPI_Group fine_grid_group;
  MPI_Group group_world;

  MPI_Comm coarse_grid_comm;
  MPI_Comm fine_grid_comm;

  int coarse_grid_ranks[] = {0,1,2,3};
  int fine_grid_ranks[] = {4};

  MPI_Comm_group(MPI_COMM_WORLD, &group_world);

  MPI_Group_incl(group_world, 4, coarse_grid_ranks, &coarse_grid_group); 
  MPI_Comm_create(MPI_COMM_WORLD, coarse_grid_group, &coarse_grid_comm);

  MPI_Group_incl(group_world, 1, fine_grid_ranks, &fine_grid_group);
  MPI_Comm_create(MPI_COMM_WORLD, fine_grid_group, &fine_grid_comm);

  int in_refine = !(rank < P*Q);


  MPI_Type_vector(par_rows - 2, 1, par_cols, MPI_DOUBLE, &MPI_CLM);
  MPI_Type_commit(&MPI_CLM);
  MPI_Type_vector(par_ref_rows - 2, 1, par_ref_cols, MPI_DOUBLE, &MPI_CLM_FINE);
  MPI_Type_commit(&MPI_CLM_FINE);

  if(!in_refine) {
    int prds[] = {0,0};
    int dims[] = {Q,P};
    MPI_Cart_create(coarse_grid_comm, 2, dims, prds, 1, &MPI_COMM_CART);
    MPI_Comm_rank(MPI_COMM_CART, &cRank);
    MPI_Cart_coords(MPI_COMM_CART, cRank, 2, crds);
  } else {
    //Create cart
    int rprds[] = {0,0};
    int rdims[] = {R,S};
    MPI_Cart_create(fine_grid_comm, 2, rdims, rprds, 1, &MPI_COMM_CART2);
    MPI_Comm_rank(MPI_COMM_CART2, &cRank2);
    MPI_Cart_coords(MPI_COMM_CART2, cRank2, 2, rcrds);
  }

  
  //width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  dt = .001;
  //dt = 5;
  int crse_rows_cutoff = refrows/2 + 1;
  int crse_cols_cutoff = refcols/2 + 1;

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
      initGrid(mat1, par_rows, par_cols, crse_rows_cutoff, crse_cols_cutoff);
      initGrid(mat2, par_rows, par_cols, crse_rows_cutoff, crse_cols_cutoff);
      //printGrid(mat1, par_rows, par_cols);
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
      //printGrid(mat1, par_ref_rows, par_ref_cols);
    }

  //arrays to store interfaces
  double *toparray, *rightarray;
  int szarray = szofmesh - 1;
  toparray = calloc(szarray, sizeof(double));
  rightarray = calloc(szarray, sizeof(double));

  //buffers to send data
  double *bfr_right_refine_pts, *bfr_top_refine_pts;
  bfr_right_refine_pts = calloc(par_ref_cols, sizeof(double));
  bfr_top_refine_pts = calloc(par_ref_rows, sizeof(double));

  time = 0;

  //Corner Points
  //coarse
  int ci_crnr = nrows - crse_rows_cutoff;
  int cj_crnr = crse_cols_cutoff - 1;
  //fine
  int fi_crnr = 0;
  int fj_crnr = refcols - 1;

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
  
  //printGrid(mat1,nrows,ncols);
  //fp = fopen ( output_file, "w" );
  //start = clock();
  for(iter = 0; iter < max_iter; iter++)
    {
      time = time + dt;

      if(!in_refine){

	computeTimestep(mat1, mat2, par_rows, par_cols, &converge, wx, wy, crse_rows_cutoff, crse_cols_cutoff);
	update(&mat1, &mat2);
      }
      else
      {
	    computeFineTimestep(mat1, mat2, refrows, refcols, &converge, wxy, sxy);
            //update(&mat1, &mat2);
      }

      if(in_refine || rank == 3)
	sendRightRefinePoints(mat2, bfr_right_refine_pts, par_ref_rows, par_ref_cols);

      if(in_refine || rank == 0)
	sendTopRefinePoints(mat2, bfr_top_refine_pts, par_ref_rows, par_ref_cols);
//recvRefinePoints
	    //Compute Interfaces
/*	    computeInterfaceRightTimestep(rightarray, mat1_refine, mat1, nrows, ncols, crse_cols_cutoff, fi_crnr, wxy, dx, crse_rows_cutoff, refcols);
	    computeInterfaceTopTimestep(toparray, mat1_refine, mat1, nrows, ncols, crse_rows_cutoff, fj_crnr, wxy, dx, crse_cols_cutoff);
	    double corner = computeCornerTimestep(mat1_refine, mat1, ci_crnr, cj_crnr, fi_crnr, fj_crnr, wxy);
	    injectCoarsePoints(toparray, rightarray, mat2, nrows, ncols, crse_rows_cutoff, crse_cols_cutoff, corner);
	    injectFinePoints(toparray, rightarray, mat2_refine, refrows, refcols, corner);
      */

      //      update(&mat1_refine, &mat2_refine, refrows, refcols);
      //if (iter%100 == 0)
      iter++;
    }

  if(!in_refine)
    printGrid(mat2, par_rows, par_cols);
  else
  {
    printGrid(mat2, par_ref_rows, par_ref_cols);
  }
 
  if (rank == 3)
    {
      for(i = 1; i < par_ref_cols - 2; i++)
	printf("right: %lf\n", bfr_right_refine_pts[i]);
    }

  if (rank == 0)
    {
      for(i = 0; i < par_ref_rows; i++)
	printf("top: %lf\n", bfr_top_refine_pts[i]);
    }
  //end = clock();
  //for(i = 0; i < szofmesh -1; i++)
  //printf("mat1 %lf\n", mat1[nrows-crse_rows_cutoff][3-1]);

  //printf("\n");
  //printGrid(mat1,nrows,ncols);
  //fprintf ( fp, "%d\n", nrows );
  //fprintf ( fp, "%d\n", ncols );


  //fprintf(fp, "num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);
  //fprintf(fp, "Total time: %f seconds\n", (((double) end ) - start)/CLOCKS_PER_SEC);
  //fclose(fp);
  //printf("num of iterations: %d\n", iter);

  MPI_Finalize();

  return 0;
}


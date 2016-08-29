
#include "grid_utils.h"
#include <stdio.h>
#include <mpi.h>

extern int P, Q, rank, size;
extern int rcrds[2], crds[2]; 

void initFineGrid(double** mat, int nrows, int ncols)
{
  int i;
  int bottom = 40;
  int left = 20;
  // left rows initial conditions
  if(rcrds[1] == 0 ) 
    {
      if(rcrds[0]%2 == 0 || nrows%2 == 0)
	{
	  for ( i = 1; i < nrows - 1; i+=2)
	    {
	      mat[i][2] = left;
	    }
	  for (i = 2; i < nrows - 1; i+=2)
	    {
	      mat[i][2] = left/2;
	    }
	}
      else 
	{
	  for ( i = 1; i < nrows - 1; i+=2)
	    {
	      mat[i][2] = left/2;
	    }
	  for (i = 2; i < nrows - 1; i+=2)
	    {
	      mat[i][2] = left;
	    }
	}
    }

  // bottom rows initial conditions
  if(rcrds[0] == Q-1)
    {
      if(rcrds[1]%2 != 1 || ncols%2 == 0) {
	for ( i = 1; i < ncols - 1; i+=2)
	  {
	    mat[nrows-3][i] = bottom/2;
	  }
	for (i = 2; i < ncols - 1; i+=2)
	  {
	    mat[nrows-3][i] = bottom;
	  }
      }
      else
	{
	  for ( i = 1; i < ncols - 1; i+=2)
	    {
	      mat[nrows-3][i] = bottom;
	    }
	  for (i = 2; i < ncols - 1; i+=2)
	    {
	      mat[nrows-3][i] = bottom/2;
	    }
	}
    } 

  // bottom-left corner is a special case - erase extra values that were put in
  if(rcrds[0]==Q-1 && rcrds[1]==0) {
    mat[nrows-3][1] = 0.0;
    mat[nrows-2][2] = 0.0;
  }
}


void initGrid(double** mat, int nrows, int ncols)
{
  int i;
  //if (rank!= 2){
  //left and right rows initial conditions
  if(crds[1] == 0)
    {
      for ( i = 1; i < nrows - 1; i++)
	{
	  mat[i][1] = 20;
	}
    }
  else if(crds[1] == P - 1)
    {
      for (i = 1; i <nrows - 1; i++)
	{ 
	  mat[i][ncols - 2] = 30;
	}
    }

  //bottom and top rows initial conditions
  if(crds[0] == 0)
    {
      for ( i = 1; i < ncols - 1; i++)
	{
	  mat[1][i] = 50;
	}
    }
  else if(crds[0] == Q-1)
    {
      for ( i = 1; i < ncols - 1; i++)
	{
	  mat[nrows-2][i] = 40;
	}
    }

  //}
}

void printGrid(double** mat, int nrows, int ncols)
{
  int i, j, k;

  for(k=0; k<size; k++) {
    if(rank==k) {
      printf("\nRank %d:\n", rank);
      for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
          printf("%.15lf ", mat[i][j]);
        }
        printf("\n");
      }
    }
  }

}


#include "ser_grid_utils.h"
#include <stdio.h>

void initFineGrid(double** mat, int nrows, int ncols)
{
  int i;
  int left = 20;
  int bottom = 40;
  //bottom and top rows initial conditions
  for ( i = 0; i < nrows; i+=2)
    {
      mat[i][0] = left;
    }
  for (i = 1; i < nrows; i+=2)
    {
      mat[i][0] = left/2;
    }
  //left and right rows initial conditions
  for ( i = 0; i < ncols; i+=2)
    {
      mat[nrows-1][i] = bottom;
    }
  for (i=1; i < ncols; i+=2)
    {
      mat[nrows-1][i] = bottom/2;
    }

}

void initGrid(double** mat, int nrows, int ncols, int crse_rows_cutoff, int crse_clms_cutoff)
{ 
  int i;
  //bottom and top rows initial conditions
  for ( i = 0; i < nrows - crse_rows_cutoff +1; i++)
    {
      mat[i][0] = 20;
    }
  for (i = 0; i <nrows; i++)
      mat[i][ncols - 1] = 30;
  //left and right rows initial conditions
  for ( i = crse_clms_cutoff -1; i < ncols; i++)
    {
      mat[nrows - 1][i] = 40;
    }
  for ( i = 0; i < ncols; i++)
      mat[0][i] = 50;

}

void printGrid(double** mat, int nrows, int ncols)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            printf("%.15lf ", mat[i][j]);
        }
        printf("\n");
    }


}




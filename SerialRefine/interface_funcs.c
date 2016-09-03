#include "interface_funcs.h"
#include <stdio.h>
void injectCoarsePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, int crse_rows_cutoff, int crse_clms_cutoff, double corner)
{
  int i, j = 0;
  //top interface row
  for (i = 1; i < crse_clms_cutoff; i++)
      mat[nrows-crse_rows_cutoff][i] = toparray[i-1];
  //right interface row
  for(i = nrows-2; i > nrows - crse_rows_cutoff; i--)
  {
      mat[i][crse_clms_cutoff - 1] = rightarray[j];
      j++;
  }
  //corner point
  mat[nrows-crse_rows_cutoff][crse_clms_cutoff-1] = corner;
}

void injectFinePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, double corner)
{
  int i, j = 0, n = 0;
  //top
  for (i = 2; i < ncols - 1; i+=2)
  {
      mat[0][i] = toparray[j];
      j++;
  }
  //right
for(i = nrows - 3; i > 0; i-=2)
  {
      mat[i][ncols - 1] = rightarray[n];
      n++;
  }
  //corner
  mat[0][ncols-1] = corner;
}

/****    Interface Points *****/
//calculates the corner point
double computeCornerTimestep(double** mat1_refine, double** mat1, int i, int j, int m, int n, double wxy)
{
  double mat2_corner;

  mat2_corner = /*prev pointmat1_corner*/mat1[i][j] + wxy * (16/15) * (/*prev point*/-4* mat1[i][j]/*mat1_corner*/+
  .5 * /*bound lwr*/mat1[i+1][j] + /*coarse*/mat1[i][j+1] + /*bound*/ mat1[i-1][j] + .5 * /*coarse*/mat1[i][j-1]
  +/*mesh*/ mat1_refine[m+1][n-1]);
 return mat2_corner;
}


void computeInterfaceRightTimestep(double* rightarray, double** mat1_refine, double** mat1, int nrows, int ncols, int crse_clms_cutoff, int m, double wxy, double dx, int crse_rows_cutoff, int refcols, int refrows)
{
  int j = crse_clms_cutoff-1, n = refcols-1, t = 0, i;
  m = refrows - 3;

    for (i = nrows-2; i > nrows-crse_rows_cutoff ; i--)
    {
      rightarray[t] =/*previous interface point*/ mat1[i][j] + wxy * (mat1[i+1][j] -2 * mat1[i][j] + mat1[i-1][j]) + wxy *
      (4/3) * (/*previous interface point*/- 4*mat1[i][j] + /*coarse*/mat1[i][j+1]  + /*coarse*/ .5 * mat1[i+1][j]\
      +/*coarse*/ .5 * mat1[i-1][j] +/*refine*/ .5 * mat1_refine[m-1][n-1] +/*refine*/ mat1_refine[m][n-1] +/*refine*/ .5\
      * mat1_refine[m+1][n-1]);//changed one n to - different from the book
      m-=2;
      //iterates through right array
      t++;
    }
}

void computeInterfaceTopTimestep(double* toparray, double** mat1_refine, double** mat1, int nrows, int ncols, int crse_rows_cutoff, int n, double wxy, double dx, int crse_clms_cutoff)
{
  int i = nrows - crse_rows_cutoff;
  int m = 0, t = 0, j;
  n = 2;


  for (j = 1; j < crse_clms_cutoff-1; j++)
  {
      toparray[t] = /*previous interface point*/mat1[i][j] + wxy * (mat1[i][j-1] -2 * mat1[i][j] + mat1[i][j+1]) + wxy * (4/3) *
      (-/*prev interface point*/ 4*mat1[i][j] +/*coarse*/ mat1[i-1][j]  +/*coarse*/ .5 * mat1[i][j+1]
      +/*coarse*/ .5 * mat1[i][j-1] +/*refine*/ .5 * mat1_refine[m+1][n-1] +/*refine*/ mat1_refine[m+1][n]
      +/*refine*/ .5 * mat1_refine[m+1][n+1]);
      n+=2;
      //iterates through top array
      t++;
  }


}


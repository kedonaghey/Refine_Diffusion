#include "compute_timesteps.h"


void computeFineTimestep(double** mat1_refine, double** mat2_refine, int nrows, int ncols, double wxy, double sxy)
{

  int i,j;

  for (i = 1; i < nrows -1; i++)
  {
    for (j = 1; j < ncols -1; j++)
    {
      //calculates cases on boundaries by interface and divides two closest boundary points
      if (i==1 && j == ncols-2)
        {
         mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i-1][j+1] \
         + mat1_refine[i-1][j-1]) -2 * mat1_refine[i][j] + mat1_refine[i+1][j]) + wxy * ( .5 * (mat1_refine[i+1][j+1]
         + mat1_refine[i-1][j+1]) -2 * mat1_refine[i][j] + mat1_refine[i][j-1]);

        }
      else if (i == 1 && j%2 == 1)
        {
          mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i-1][j+1] \
         + mat1_refine[i-1][j-1]) -2 * mat1_refine[i][j] + mat1_refine[i+1][j]) + wxy * (mat1_refine[i][j-1] -2 * mat1_refine[i][j] \
         + mat1_refine[i][j+1]);
        }
      else if (j == ncols-2 && i%2 == 1)
        {
         mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i+1][j+1]
         + mat1_refine[i-1][j+1]) -2 * mat1_refine[i][j] + mat1_refine[i][j-1]) + wxy * (mat1_refine[i+1][j] -2 * mat1_refine[i][j]\
         + mat1_refine[i-1][j]);
        }
      else
        {
         mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]\
         -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
        }
    }
  }

}

void computeTimestep(double** mat1, double** mat2, int nrows, int ncols, double wx, double wy)
{
  int i, j;
 
  for (i = 1; i < nrows - 1; i++)
  {
    for (j = 1; j < ncols - 1; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
    }
  }
}








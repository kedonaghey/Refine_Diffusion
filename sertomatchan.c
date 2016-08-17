#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void InitGrid(double** mat, int nrows, int ncols)
{
  int i, j;
  //bottom and top rows initial conditions
  for ( i = 0; i < nrows; i++)
    {
      mat[i][0] = 0;//20;//0
      mat[i][ncols - 1] = 0;//30;//100
    }
  //left and right rows initial conditions
  for ( i = 0; i < ncols; i++)
    {
      mat[0][i] = 0; //50
      mat[nrows - 1][i] = 0;//40;//0
    }

  //
  for(i = 1; i < nrows -1; i++)
  {
    for(j = 1; j < ncols -1; j++)
    {
      mat[i][j] = 50;
    }
  }
/*  for(i = ncols/2; i < nrows - 1; i++)
  {
    for(j = 1; j < ncols - 1; j++)
    {
      mat[i][j] = 50;
    }
  }
*/
}

void Printgrid(double** mat, int nrows, int ncols)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }


}
void PrintGrid(double** mat, int nrows, int ncols, FILE *fp)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            fprintf(fp, "%f ", mat[i][j]);
        }
        fprintf(fp, "\n");
    }


}
void Compute_TimeStep(double** mat1, double** mat2, int nrows, int ncols, double* converge, double dt, double dx, double dy)
{
  int i, j;
  double **tmp;
  double wx = 0, wy = 0;
  //diffusitivity coefficient
  double k = 1;

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

  
  double diff = 0;
  for (i = 1; i < nrows - 1; i++)
  {
    for (j = 1; j < ncols - 1; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( diff < fabs(mat2[i][j] - mat1[i][j]))
      {
        diff = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }

  	*converge = diff;  

  for(i = 0; i < nrows; i++)
  {
    for(j = 0; j < ncols; j++)
    {
      mat1[i][j] = mat2[i][j];
    }
  }

//  Printgrid(mat1,nrows,ncols);
  //swap mat1 = mat2
/*tmp = mat1;
  mat1 = mat2;
  mat2 = tmp;*/
}


int main(int argc, char *argv[])
{

  int i, j;
  double **mat1, **mat2;
  //size of grid
  int nrows = 11, ncols = 11;
  double dx = 0, dy= 0, dt, time;
  double converge = 0, epsilon;
  int iter;
  int max_iter = atoi(argv[1]);
  clock_t start, end;
  int prnt = 1;
  char output_file[80] = "grid.dat";
  FILE *fp;

  //width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  epsilon = .001;
  dt = .0001;
  //dt = 5;
  mat1 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat1[i] = calloc(ncols, sizeof(double));
	

  mat2 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat2[i] = calloc(ncols, sizeof(double));

  InitGrid(mat1, nrows, ncols);
  InitGrid(mat2, nrows, ncols);

  time = 0;

  fp = fopen ( output_file, "w" );
  start = clock();
  for(iter = 0; iter < max_iter; iter++)
  {
    //move to seperate func?
    //while ( converge >= epsilon)
    //{
      time = time + dt;
      Compute_TimeStep(mat1, mat2, nrows, ncols, &converge, dt, dx, dy);
      iter++;
      if (converge < epsilon)
       	break;
    //}
  }
  end = clock();

 Printgrid(mat1,nrows,ncols);
  fprintf ( fp, "%d\n", nrows );
  fprintf ( fp, "%d\n", ncols );
/*
  for(j=ncols+1;j>=0;j--)
  { for (i=0;i<=nrows;i++)
    { fprintf(fp,"%15.11f ",x0[j][i]);}
    fprintf(fp,"%15.11f\n",x0[j][nrows+1]);
  }
*/
  fprintf(fp, "num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);
  fprintf(fp, "Total time: %f seconds\n", (((double) end ) - start)/CLOCKS_PER_SEC);

  fclose(fp);
  printf("num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);
}


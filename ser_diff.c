#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

void initGrid(double** mat, int nrows, int ncols)
{
  int i;
  //bottom and top rows initial conditions
  for ( i = 0; i < nrows; i++)
    {
      mat[i][0] = 20;//0
      mat[i][ncols - 1] = 30;//100
    }
  //left and right rows initial conditions
  for ( i = 0; i < ncols; i++)
    {
      mat[0][i] = 50;//100
      mat[nrows - 1][i] = 40;//0
    }

}

void printGrid(double** mat, int nrows, int ncols)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }


}

void printGridToFile(double** mat, int nrows, int ncols, FILE *fp)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            fprintf(fp, "%f ", mat[i][j]);
        }
        fprintf(fp, "\n");
    }


}
void computeTimestep(double*** mat1_ptr, double*** mat2_ptr, int nrows, int ncols, double* converge, double wx, double wy)
{
  int i, j;
  double **mat1 = *mat1_ptr;
  double **mat2 = *mat2_ptr;

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
 
  //swap mat1 = mat2 
  double** tmp = *mat1_ptr;
  *mat1_ptr = *mat2_ptr;
  *mat2_ptr = tmp;
}


int main(int argc, char *argv[])
{

  int i, c;
  double **mat1, **mat2;
  //size of grid
  int nrows = 16, ncols = 16;
  double dx = 0, dy= 0, dt, time;
  double converge = 0, epsilon;
  int iter, max_iter = 10;
  //clock_t start, end;
  int prnt = 1;
  char output_file[80] = "heat.dat";
  FILE *fp;

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


  //width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  epsilon = .001;
  dt = .00001;
  //dt = 5;
  mat1 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat1[i] = calloc(ncols, sizeof(double));
	

  mat2 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat2[i] = calloc(ncols, sizeof(double));

  initGrid(mat1, nrows, ncols);
  initGrid(mat2, nrows, ncols);

  //printGrid(mat1,nrows,ncols);
  time = 0;


  //calculate weights
  double wx = 0, wy = 0;
  //diffusitivity coefficient
  double k = 1;

//printf("dx = %lf, dy = %lf, dt = %lf\n", dx, dy, dt);
  //weight in x direction
  wx = k * (dt/(dx*dx));
  //weight in y direction
  wy = k * (dt/(dy*dy));

  //stability
  // wx + wy must be less than .5
  if ( wx + wy > .5)
  {
    printf("Does not reach requirements for stability\n");
    printf("wx = %lf, wy = %lf\n", wx, wy);
    exit(1);
  }





  fp = fopen ( output_file, "w" );
  //start = clock();
  for(iter = 0; iter < max_iter; iter++)
  {
    //move to seperate func?
    //while ( converge >= epsilon)
    //{
      time = time + dt;
      computeTimestep(&mat1, &mat2, nrows, ncols, &converge, wx, wy);
      if (iter%100 == 0)
      { 
        //printf("%d %f\n", iter, converge);
	    //prnt = 2 * prnt;
      }
     
      if (converge < epsilon)
       	break;
    //}

      printf("%lf\n", converge);
  }
  //end = clock();

  //fprintf ( fp, "%d\n", nrows );
  //fprintf ( fp, "%d\n", ncols );

  //printGridToFile(mat1, nrows, ncols, fp);
  //fprintf(fp, "num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);
  //fprintf(fp, "Total time: %f seconds\n", (((double) end ) - start)/CLOCKS_PER_SEC);

  fclose(fp);
  printf("num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);

  return 0;
}


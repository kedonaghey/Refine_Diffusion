#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void InitFineGrid(double** mat, int nrows, int ncols)
{
  int i;
  int top = 20;
  int left = 50;
  //bottom and top rows initial conditions
  for ( i = 0; i < nrows; i+=2)
    {
      mat[i][0] = top;
    }
  for (i = 1; i < nrows; i+=2)
    {
      mat[i][0] = top/2;
    }
  //left and right rows initial conditions
  for ( i = 0; i < ncols; i+=2)
    {
      mat[nrows-1][i] = left;//100
    }
  for (i=1; i < ncols; i+=2)
    {
      mat[nrows-1][i] = left/2;
    }

}

void InitGrid(double** mat, int nrows, int ncols, int crs, int ccs)
{
  int i;
  //bottom and top rows initial conditions
  for ( i = 0; i < nrows - crs +1; i++)
    {
      mat[i][0] = 20;//0
    }
  for (i = 0; i <nrows; i++)
      mat[i][ncols - 1] = 30;//100
  //left and right rows initial conditions
  for ( i = ccs -1; i < ncols; i++)
    {
      mat[nrows - 1][i] = 50;//0
    }
  for ( i = 0; i < ncols; i++)
      mat[0][i] = 40;//100

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

void InjectCoarsePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, int crs, int ccs, double corner)
{
  int i, j;
  int test[3] = {1,2,3};
  //top
  for (i = 1; i < ccs; i++)
      mat[nrows-crs][i] = test[i-1];//toparray[i-1];//0
  //right
  for(i = nrows-2; i > nrows - crs; i--)
      mat[i][ccs - 1] = test[i-crs-1];//rightarray[i-crs-1];//100
  //corner
  mat[nrows-crs][ccs-1] = corner;
}

void InjectFinePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, double corner)
{
  int i, j;
  int test[7] = {1,2,3,4,5,6,7};
  //top
  for (i = 1; i < ncols - 1; i++)
      mat[0][i] = test[i-1];//toparray[i-1];//0
  //right
  for(i = 1; i < nrows - 1; i++)
      mat[i][ncols - 1] = test[i-1];//rightarray[i-1];//100
  //corner
  //double check
  mat[0][ncols-1] = corner;
}

/****    Interface Points *****/
//for bound points use regular matrix

double Compute_CornerTimeStep(double** mat1_refine, double** mat1, int i, int j, int m, int n, double dt, double dx, double dy)
{
  double mat1_corner;
  double mat2_corner;
  double wxy = 0;
  double k = .02, diff = 0;
  wxy = k * dt/(dx*dx);

  mat2_corner = /*prev point*/mat1_corner + wxy * (16/15) * (/*prev point*/-4*mat1_corner+ .5 * /*bound*/mat1[i][j-1]\
  + /*coarse*/mat1[i+1][j] + .5 * /*coarse*/mat1[i-1][j] +/*bound*/ mat1[i][j+1] +/*mesh*/ mat1_refine[m+1][n-1]);

 return mat2_corner;
}

//need to add in different

void Compute_InterfaceRightTimeStep(double* rightarray, double** mat1_refine, double** mat1, int nrows, int ncols, int i, int m, double dt, double dx, double dy, int crs)
{
  int j, n = 1, t = 0;
  double wxy = 0;
  double k = .02, diff = 0;
  wxy = k * dt/(dx*dx);

  //actually probably dont need
  //if (i,j-1) from corner or (i-1,j) from corner - substitute in corner
  //else do the following
  //do one for loop hardcode i,m stays same
    for (j = nrows-2; j > nrows-crs; j--)
    {
      rightarray[t] =/*previous interface point*/ mat1[i][j] + /*previous interface point*/ mat1[i][j] * wxy * (dx/2) + wxy * (4/3) * \
      (/*previous interface point*/- 4*mat1[i][j] + /*coarse*/mat1[i+1][j]  + /*coarse*/ .5 * mat1[i][j+1]\
      +/*coarse*/ .5 * mat1[i][j-1] +/*refine*/ .5 * mat1_refine[m+1][n+1] +/*refine*/ mat1_refine[m+1][n] +/*refine*/ .5 * mat1_refine[m+1][n-1]);//changed one n to - dif\
      ferent from the book
      n+=2;
      t++;
    }


}


void Compute_InterfaceTopTimeStep(double* toparray, double** mat1_refine, double** mat1, int nrows, int ncols, int j, int n,double dt, double dx, double dy, int ccs)
{
  int i, m = 1, t = 0;
  double wxy = 0;
  double k = .02, diff = 0;
  wxy = k * dt/(dx*dx);


  //actually probably dont need
  //if (i,j-1) from corner or (i-1,j) from corner - substitute in corner
  //else do the following
  //do one for loop hardcode j,n stays same
  //after declare m/n += 2 ??
  for (i = 1; i < ccs; i++)
  {
      toparray[t] = /*previous interface point*/mat1[i][j] + /*prev interface point*/mat1[i][j] * wxy * (dx/2) + wxy * (4/3) * \
      (-/*prev interface point*/ 4*mat1_refine[i][j] +/*coarse*/ mat1[i][j+1]  +/*coarse*/ .5 * mat1_refine[i+1][j]\
      +/*coarse*/ .5 * mat1_refine[i-1][j] +/*refine*/ .5 * mat1_refine[m+1][n-1] +/*refine*/ mat1_refine[m][n-1] +/*refine*/ .5 * mat1_refine[m-1][n-1]);
      m+=2;
      t++;
  }


}


void Compute_FineTimeStep(double** mat1_refine, double** mat2_refine, int nrows, int ncols, double* converge, double dt, double dx, double dy)
{

  int i,j;
  double wxy = 0;
  double k = .02, diff = 0;
  int refinement = 2;
  double sxy = (double)dx/refinement;
  wxy = k * dt/(dx*dx);

  for (i = 1; i < nrows -1; i++)
  {
    for (j = 1; j < ncols -1; j++)
    {
      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i+1][j+1] + mat1_refine[i+1][j-1])\
      -2 * mat1_refine[i][j] + mat1_refine[i-1][j]);
/*      if ( diff < fabs(mat2_refine[i][j] - mat1_refine[i][j]))
      {
        diff = fabs(mat1_refine[i][j] - mat2_refine[i][j]);
      }
    }
  }

  *converge = diff;  
*/
}}
/*  for(i = 0; i < nrows; i++)
  {
    for(j = 0; j < ncols; j++)
    {
      mat1_refine[i][j] = mat2_refine[i][j];
    }
  }
*/
}

void Compute_TimeStep(double** mat1, double** mat2, int nrows, int ncols, double* converge, double dt, double dx, double dy)
{
  int i, j;
  double **tmp;
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

  
  double diff = 0;
  for (i = 1; i < nrows -1; i++)
  {
    for (j = 1; j < ncols -1; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
/*      if ( diff < fabs(mat2[i][j] - mat1[i][j]))
      {
        diff = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }

  	*converge = diff;  
*/}}
/*  for(i = 0; i < nrows; i++)
  {
    for(j = 0; j < ncols; j++)
    {
      mat1[i][j] = mat2[i][j];
    }
  }
*/
}

void Update(double** mat1, double** mat2, int nrows, int ncols)
{
  int i, j;
  for(i = 0; i < nrows; i++)
  {
    for(j = 0; j < ncols; j++)
      mat1[i][j] = mat2[i][j];
  }
}

int main(int argc, char *argv[])
{

  int i, j;
  double **mat1, **mat2;
  double **mat1_refine, **mat2_refine;
  //size of grid
  int nrows = 10, ncols = 10;
  //2x2 for mesh in a 4x4 array
  int szofmesh = 4;
  int refrows = szofmesh * 2 + 1;
  int refcols = szofmesh * 2 + 1;
  double dx = 0, dy= 0, dt, time;
  double converge = 0, epsilon;
  int iter, max_iter = 100;
  clock_t start, end;
  int prnt = 1;
  char output_file[80] = "k.txt";
  FILE *fp;

  //width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  epsilon = .001;
  dt = .001;
  //dt = 5;
  mat1 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat1[i] = calloc(ncols, sizeof(double));
	

  mat2 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat2[i] = calloc(ncols, sizeof(double));

  int crs = refrows/2 + 1;
  int ccs = refcols/2 + 1;
  InitGrid(mat1, nrows, ncols, crs, ccs);
  InitGrid(mat2, nrows, ncols, crs, ccs);

  mat1_refine = calloc(refrows, sizeof(double *));
  for (i = 0; i < refrows; i++)
    mat1_refine[i] = calloc(refcols, sizeof(double));
	

  mat2_refine = calloc(refrows, sizeof(double *));
  for (i = 0; i < refrows; i++)
    mat2_refine[i] = calloc(refcols, sizeof(double));

  InitFineGrid(mat1_refine, refrows, refcols);
  InitFineGrid(mat2_refine, refrows, refcols); 

  //arrays to store interfaces
  double *toparray, *rightarray;
  int szarray = szofmesh - 1;
  toparray = calloc(szarray, sizeof(double));
  rightarray = calloc(szarray, sizeof(double));

  time = 0;

  //Corner Points
  //coarse
  int ci = nrows - crs;
  int cj = ccs - 1;
  //fine
  int fi = 0;
  int fj = refcols - 1;

  //Printgrid(mat1,nrows,ncols);
  fp = fopen ( output_file, "w" );
  start = clock();
  for(iter = 0; iter < max_iter; iter++)
  {
    //move to seperate func?
    //while ( converge >= epsilon)
    //{
      time = time + dt;
      Compute_TimeStep(mat1, mat2, nrows, ncols, &converge, dt, dx, dy);
      Compute_FineTimeStep(mat1_refine, mat2_refine, refrows, refcols, &converge, dt, dx, dy);
      //Compute Interfaces
      Compute_InterfaceRightTimeStep(rightarray, mat1_refine, mat1, nrows, ncols, ci, fi, dt, dx, dy, crs);
      Compute_InterfaceTopTimeStep(toparray, mat1_refine, mat1, nrows, ncols, cj, fj, dt, dx, dy, ccs);
      double corner = Compute_CornerTimeStep(mat1_refine, mat1, ci, cj, fi, fj, dt, dx, dy);
      InjectCoarsePoints(toparray, rightarray, mat2, nrows, ncols, crs, ccs, 100/*corner*/);
      InjectFinePoints(toparray, rightarray, mat2_refine, refrows, refcols, 100/*corner*/);
      Update(mat1, mat2, nrows, ncols);
      Update(mat1_refine, mat2_refine, refrows, refcols);
      //if (iter%100 == 0)
      //{ 
      //  printf("%d %f\n", iter, converge);
      //  fprintf(fp, "\n");
	    //prnt = 2 * prnt;
      //}
      iter++;
  //    if (converge < epsilon)
    //   	break;
    //}
  }
  end = clock();

  mat1[0][0] = 123456;
  Printgrid(mat1_refine,refrows,refcols);
  printf("\n");
  Printgrid(mat1,nrows,ncols);
  //fprintf ( fp, "%d\n", nrows );
  //fprintf ( fp, "%d\n", ncols );

  PrintGrid(mat1, nrows, ncols, fp);
  //fprintf(fp, "num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);
  //fprintf(fp, "Total time: %f seconds\n", (((double) end ) - start)/CLOCKS_PER_SEC);

  fclose(fp);
  printf("num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);


}


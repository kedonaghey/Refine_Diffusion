#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

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
      mat[nrows-1][i] = bottom;//100
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
      mat[i][0] = 20;//0
    }
  for (i = 0; i <nrows; i++)
      mat[i][ncols - 1] = 30;//100
  //left and right rows initial conditions
  for ( i = crse_clms_cutoff -1; i < ncols; i++)
    {
      mat[nrows - 1][i] = 40;//0
    }
  for ( i = 0; i < ncols; i++)
      mat[0][i] = 50;//100

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
void printGridToFile(double** mat, int nrows, int ncols, FILE *fp)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            fprintf(fp, "%8f ", mat[i][j]);
        }
        fprintf(fp, "\n");
    }


}

void injectCoarsePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, int crse_rows_cutoff, int crse_clms_cutoff, double corner)
{
  int i, j = 0;
  //int test[] = {1,2,3};
  //top
  for (i = 1; i < crse_clms_cutoff; i++)
      mat[nrows-crse_rows_cutoff][i] = toparray[i-1];//toparray[i-1];//0
  //right
  for(i = nrows-2; i > nrows - crse_rows_cutoff; i--)
  {
      //mat[i][crse_clms_cutoff - 1] = test[i-crse_rows_cutoff-1];//rightarray[i-crs-1];//100
      mat[i][crse_clms_cutoff - 1] = rightarray[j];
      j++;
  }
  //corner
  mat[nrows-crse_rows_cutoff][crse_clms_cutoff-1] = corner;
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
  printf("mat[i][j] %.15lf\n mat1[i+1][j] %.15lf\n mat1[i][j+1] %.15lf\n mat1[i-1][j] %.15lf\n mat1[i][j-1] %.15lf\n mat1refine[i+1][j] %.15lf\n", mat1[i][j], mat1[i+1][j], mat1[i][j+1], mat1[i-1][j], mat1[i][j-1], mat1_refine[m+1][n-1]);
printf("updated value: %.15lf\n", mat2_corner);
 return mat2_corner;
}

//need to add in different

void computeInterfaceRightTimestep(double* rightarray, double** mat1_refine, double** mat1, int nrows, int ncols, int crse_clms_cutoff, int m, double wxy, double dx, int crse_rows_cutoff, int refcols, int refrows)
{
  int j = crse_clms_cutoff-1, n = refcols-1, t = 0, i;
  m = refrows - 3;

    for (i = nrows-2; i > nrows-crse_rows_cutoff ; i--)
    {
      rightarray[t] =/*previous interface point*/ mat1[i][j] + /*previous interface point*/ mat1[i][j] * wxy * (dx/2) + wxy *\
      (4/3) * (/*previous interface point*/- 4*mat1[i][j] + /*coarse*/mat1[i][j+1]  + /*coarse*/ .5 * mat1[i+1][j]\
      +/*coarse*/ .5 * mat1[i-1][j] +/*refine*/ .5 * mat1_refine[m-1][n-1] +/*refine*/ mat1_refine[m][n-1] +/*refine*/ .5\
      * mat1_refine[m+1][n-1]);//changed one n to - different from the book
      m-=2;
      //iterates through right array
      t++;
    }

    for(i=0; i<4; i++)
      printf("%lf\n", rightarray[i]);
    printf("\n");


}


void computeInterfaceTopTimestep(double* toparray, double** mat1_refine, double** mat1, int nrows, int ncols, int crse_rows_cutoff, int n, double wxy, double dx, int crse_clms_cutoff)
{
  int i = nrows - crse_rows_cutoff;
  int m = 0, t = 0, j;
  n = 2;


  for (j = 1; j < crse_clms_cutoff-1; j++)
  {
      printf("wxy %lf, dx/2 %lf\n", wxy, dx/2);
      toparray[t] = /*previous interface point*/mat1[i][j] + /*prev interface point*/mat1[i][j] * wxy * (dx/2) + wxy * (4/3) *
      (-/*prev interface point*/ 4*mat1[i][j] +/*coarse*/ mat1[i-1][j]  +/*coarse*/ .5 * mat1[i][j+1]
      +/*coarse*/ .5 * mat1[i][j-1] +/*refine*/ .5 * mat1_refine[m+1][n-1] +/*refine*/ mat1_refine[m+1][n] 
      +/*refine*/ .5 * mat1_refine[m+1][n+1]);
//      printf("i = %d, j = %d, m = %d n = %d\n", i, j, m , n);
//      printf("mat[i][j] %.15lf\n mat1[i-1][j] %.15lf\n mat1[i][j+1] %.15lf\n mat1[i][j-1] %.15lf\n mat1[i+1][j-1] %.15lf\n mat1[i+1][j] %.15lf\n mat[i+1][j+1] %.15lf\n", mat1[i][j], mat1[i-1][j], mat1[i][j+1], mat1[i][j-1], mat1_refine[m+1][n-1], mat1_refine[m+1][n], mat1_refine[m+1][n+1]);
//      printf("updated value: %.15lf\n", toparray[t]);
      n+=2;
      //iterates through top array
      t++;
  }

//  printf("%lf\n", toparray[1]);
//  printf("%lf\n", toparray[2]);

}


void computeFineTimestep(double** mat1_refine, double** mat2_refine, int nrows, int ncols, double* converge, double wxy, double sxy)
{

  int i,j;

  for (i = 1; i < nrows -1; i++)
  {
    for (j = 1; j < ncols -1; j++)
    {
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


void computeTimestep(double** mat1, double** mat2, int nrows, int ncols, double* converge, double wx, double wy, int crse_rows_cutoff, int crse_clms_cutoff)
{
  int i, j;
 
  for (i = 1; i < nrows - 1; i++)
  {
    for (j = 1; j < ncols - 1; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
}}

}

void update(double*** mat1_ptr, double*** mat2_ptr, int nrows, int ncols)
{
  double **tmp = *mat1_ptr;
  *mat1_ptr = *mat2_ptr;
  *mat2_ptr = tmp;
}

int main(int argc, char *argv[])
{

  int i, c;
  double **mat1, **mat2;
  double **mat1_refine, **mat2_refine;
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

  int szofmesh = .5*nrows;

  int refrows = szofmesh * 2 - 1;
  int refcols = szofmesh * 2 - 1;
  //width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  dt = .001;
  //dt = 5;
  mat1 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat1[i] = calloc(ncols, sizeof(double));
	

  mat2 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat2[i] = calloc(ncols, sizeof(double));

  int crse_rows_cutoff = ceil(refrows/2.0) + 0;
  int crse_clms_cutoff = ceil(refcols/2.0) + 0;
  initGrid(mat1, nrows, ncols, crse_rows_cutoff, crse_clms_cutoff);
  initGrid(mat2, nrows, ncols, crse_rows_cutoff, crse_clms_cutoff);

  mat1_refine = calloc(refrows, sizeof(double *));
  for (i = 0; i < refrows; i++)
    mat1_refine[i] = calloc(refcols, sizeof(double));
	

  mat2_refine = calloc(refrows, sizeof(double *));
  for (i = 0; i < refrows; i++)
    mat2_refine[i] = calloc(refcols, sizeof(double));

  initFineGrid(mat1_refine, refrows, refcols);
  initFineGrid(mat2_refine, refrows, refcols); 

  //arrays to store interfaces
  double *toparray, *rightarray;
  int szarray = szofmesh - 1;
  toparray = calloc(szarray, sizeof(double));
  rightarray = calloc(szarray, sizeof(double));

  time = 0;

  //Corner Points
  //coarse
  int ci_crnr = nrows - crse_rows_cutoff;
  int cj_crnr = crse_clms_cutoff - 1;
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
  fp = fopen ( output_file, "w" );
  //start = clock();
  for(iter = 0; iter < max_iter; iter++)
  {
      time = time + dt;
      computeTimestep(mat1, mat2, nrows, ncols, &converge, wx, wy, crse_rows_cutoff, crse_clms_cutoff);
      computeFineTimestep(mat1_refine, mat2_refine, refrows, refcols, &converge, wxy, sxy);
      //Compute Interfaces
      computeInterfaceRightTimestep(rightarray, mat1_refine, mat1, nrows, ncols, crse_clms_cutoff, fi_crnr, wxy, dx, crse_rows_cutoff, refcols, refrows);
      computeInterfaceTopTimestep(toparray, mat1_refine, mat1, nrows, ncols, crse_rows_cutoff, fj_crnr, wxy, dx, crse_clms_cutoff);
      double corner = computeCornerTimestep(mat1_refine, mat1, ci_crnr, cj_crnr, fi_crnr, fj_crnr, wxy);
      injectCoarsePoints(toparray, rightarray, mat2, nrows, ncols, crse_rows_cutoff, crse_clms_cutoff, corner);
      injectFinePoints(toparray, rightarray, mat2_refine, refrows, refcols, corner);
      update(&mat1, &mat2, nrows, ncols);
      update(&mat1_refine, &mat2_refine, refrows, refcols);
      //if (iter%100 == 0)
    }
  //end = clock();
  //for(i = 0; i < szofmesh -1; i++)
  //printf("mat1 %lf\n", mat1[nrows-crse_rows_cutoff][3-1]);
  printf("fine %lf\n", mat1_refine[1][refcols-2]);

  printGrid(mat1_refine,refrows,refcols);
  printf("\n");
  printGrid(mat1,nrows,ncols);
  //fprintf ( fp, "%d\n", nrows );
  //fprintf ( fp, "%d\n", ncols );

  printGridToFile(mat1, nrows, ncols, fp);
  //fprintf(fp, "num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);
  //fprintf(fp, "Total time: %f seconds\n", (((double) end ) - start)/CLOCKS_PER_SEC);

  fclose(fp);
  printf("num of iterations: %d\n", iter);

  return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include "ser_grid_utils.h"
#include "interface_funcs.h"
#include "compute_timesteps.h"

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
  clock_t start, end;
  int nrows = 5000, ncols = 5000;
  double dx = 0, dy= 0, dt, time;
  double converge = 0;
  int iter, max_iter = 100;

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
  //size of refinement
  int refrows = szofmesh * 2 - 1;
  int refcols = szofmesh * 2 - 1;
  //width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  dt = .001;
  mat1 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat1[i] = calloc(ncols, sizeof(double));
	

  mat2 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat2[i] = calloc(ncols, sizeof(double));

  //sets where refinement ends on coarse grid
  int crse_rows_cutoff = ceil(refrows/2.0) + 0;
  int crse_clms_cutoff = ceil(refcols/2.0) + 0;

  //initializes boundary points on coarse grid
  initGrid(mat1, nrows, ncols, crse_rows_cutoff, crse_clms_cutoff);
  initGrid(mat2, nrows, ncols, crse_rows_cutoff, crse_clms_cutoff);

  mat1_refine = calloc(refrows, sizeof(double *));
  for (i = 0; i < refrows; i++)
    mat1_refine[i] = calloc(refcols, sizeof(double));
	

  mat2_refine = calloc(refrows, sizeof(double *));
  for (i = 0; i < refrows; i++)
    mat2_refine[i] = calloc(refcols, sizeof(double));

  //initializes boundary points on fine grid
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
  double rxy = k * dt/(sxy*sxy);

  //start = clock();
  printf("rxy: %lf\n",rxy);
  for(iter = 0; iter < max_iter; iter++)
  {
      time = time + dt;
      computeTimestep(mat1, mat2, nrows, ncols, wx, wy);
      computeFineTimestep(mat1_refine, mat2_refine, refrows, refcols, rxy);
      //Compute Interfaces
      computeInterfaceRightTimestep(rightarray, mat1_refine, mat1, nrows, ncols, crse_clms_cutoff, fi_crnr, wxy, dx, crse_rows_cutoff, refcols, refrows);
      computeInterfaceTopTimestep(toparray, mat1_refine, mat1, nrows, ncols, crse_rows_cutoff, fj_crnr, wxy, dx, crse_clms_cutoff);
      double corner = computeCornerTimestep(mat1_refine, mat1, ci_crnr, cj_crnr, fi_crnr, fj_crnr, wxy);
      injectCoarsePoints(toparray, rightarray, mat2, nrows, ncols, crse_rows_cutoff, crse_clms_cutoff, corner);
      injectFinePoints(toparray, rightarray, mat2_refine, refrows, refcols, corner);
      update(&mat1, &mat2, nrows, ncols);
      update(&mat1_refine, &mat2_refine, refrows, refcols);
    }
  //end = clock();
  printGrid(mat1_refine,refrows,refcols);
  //printf("\n");
  printGrid(mat1,nrows,ncols);
  //printf("Total time: %f seconds\n", (((double) end ) - start)/CLOCKS_PER_SEC);
  //printf("%.15lf\n", mat1_refine[refrows-2][1]);
  //printf("num of iterations: %d\n", iter);

  return 0;
}


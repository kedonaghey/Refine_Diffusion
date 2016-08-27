#ifndef __INTERFACE_UTILS_H__
#define __INTERFACE_UTILS_H__

void injectCoarsePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, int crse_rows_cutoff, int crse_clms_cutoff, double corner);
void injectFinePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, double corner);
double computeCornerTimestep(double** mat1_refine, double** mat1, int i, int j, int m, int n, double wxy);
void computeInterfaceRightTimestep(double* rightarray, double** mat1_refine, double** mat1, int nrows, int ncols, int crse_clms_cutoff, int m, double wxy, double dx, int crse_rows_cutoff, int refcols, int refrows);
void computeInterfaceRightTimestep(double* rightarray, double** mat1_refine, double** mat1, int nrows, int ncols, int crse_clms_cutoff, int m, double wxy, double dx, int crse_rows_cutoff, int refcols, int refrows);
void computeInterfaceTopTimestep(double* toparray, double** mat1_refine, double** mat1, int nrows, int ncols, int crse_rows_cutoff, int n, double wxy, double dx, int crse_clms_cutoff);

#endif

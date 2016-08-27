
#ifndef __COMPUTE_TIMESTEP_H__
#define __COMPUTE_TIMESTEP_H__

void computeCornerTimestep(double* buffer, double** mat1, double** mat2, double wxy, int par_ref_cols);
void computeInterfaceRightTimestep(double* buffer, double** mat1, double** mat2, int par_ref_rows, int ncols, double wxy, double dx, int refcols);
void computeInterfaceTopTimestep(double* buffer, double** mat1, double** mat2, int par_rows, int par_ref_cols, double wxy, double dx);
void computeFineTimestep(double** mat1_refine, double** mat2_refine, int par_ref_rows, int par_ref_cols, double wxy, double sxy);
void computeTimestep(double* bfr_top_refine_pts, double* bfr_right_refine_pts, double** mat1, double** mat2, int par_rows, int par_cols, double wx, double wy);

#endif

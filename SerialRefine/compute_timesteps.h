#ifndef __INTERFACE_UTILS_H__
#define __INTERFACE_UTILS_H__

void computeFineTimestep(double** mat1_refine, double** mat2_refine, int nrows, int ncols, double wxy, double sxy);
void computeTimestep(double** mat1, double** mat2, int nrows, int ncols, double wx, double wy);

#endif

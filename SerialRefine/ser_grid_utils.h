#ifndef __SER_GRID_UTILS_H__
#define __SER_GRID_UTILS_H__


void initFineGrid(double** mat, int nrows, int ncols);
void initGrid(double** mat, int nrows, int ncols, int crse_rows_cutoff, int crse_clms_cutoff);
void printGrid(double** mat, int nrows, int ncols);
//void printGridToFile(double** mat, int nrows, int ncols, FILE *fp);

#endif

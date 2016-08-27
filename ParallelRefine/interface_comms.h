#ifndef __INTERFACE_COMMS_H__
#define __INTERFACE_COMMS_H__

void injectCoarsePoints(double* toparray, double* rightarray, double** mat, int par_rows, int par_cols);
void injectFinePoints(double* topbuffer, double* rightbuffer, double** mat, int par_ref_rows, int par_cols, int par_ref_cols);
void sendRightPoints(double** mat, double* buffer, int par_rows, int par_cols);
void sendTopPoints(double** mat, double* buffer, int par_rows, int par_ref_cols);
void sendRightInterfacePoints(double** mat, double* buffer, int par_ref_rows, int par_ref_cols);
void sendTopInterfacePoints(double** mat, double* buffer, int par_ref_rows, int par_ref_cols);

#endif

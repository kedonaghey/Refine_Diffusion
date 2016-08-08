#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

int rank, size;
int P, Q;
int rows, cols;

void InitGrid(double** mat)
{
  int i;
  if(rank%P == 0)
  {
    for(i= 0;i<rows;i++)
	mat[i][1] = 20;
  } 
  else if(rank%P == P-1)
  {
    for(i=0;i<rows ;i++)
	mat[i][cols-2] = 30;
  }

  if(rank/P == 0)
  {
    for(i=0 ;i<cols;i++)
	mat[1][i] = 50;
  }
  else if(rank/P == Q-1)
  {
    for(i=0;i<cols;i++)
	mat[rows - 2][i] = 40;
  }
}

void PrintGrid(double** mat, int nrows, int ncols)
{
  
  int i, j, k, l;/* 
  MPI_Barrier(MPI_COMM_WORLD);

  for(l=0; l<Q; l++) {
    if(rank/Q == l) {
      for (i = 1; i < rows-1; i++) {
	for(k=0; k<P; k++) {
	  if(rank%P == k) {
	    for (j = 1; j < cols-1; j++) {
	      printf("%f ", mat[i][j]);
	    }
	    fflush(stdout);
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	}
      }

      if(!rank) printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }  

  MPI_Barrier(MPI_COMM_WORLD);*/
/*  MPI_Barrier(MPI_COMM_WORLD);

  for(l=0; l<P; l++) {
    if(rank/P == l) {
      for (i = 1; i < rows-1; i++) {
	for(k=0; k<Q; k++) {
	  if(rank/Q == k) {
	    for (j = 1; j < cols-1; j++) {
	      printf("%f ", mat[i][j]);
	    }
            fflush(stdout);
	    printf("\n");
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	}
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }  

  MPI_Barrier(MPI_COMM_WORLD);*/
  double *line_grid;
  if (rank == 0){
     i = P*Q*rows* cols;
     line_grid = malloc(i*sizeof(double));
  }

  MPI_Gather(&mat[0][0],rows*cols,MPI_DOUBLE,line_grid,rows*cols,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (rank == 0){
  for(j=0;j<i;j++)
    printf("j =%d, %f\n ",j, line_grid[j]);
}
}

  
void Compute_TimeStep(double** mat1, double** mat2, int nrows, int ncols, double* converge, double dt, double dx, double dy, double* all_local_converge)
{
  int i, j, m;
  double wx = 0, wy = 0;
  //diffusitivity coefficient - can equal .002
  double k = 1;
  double d1 = 0.0, d2 = 0.0, d3 = 0.0, d4 = 0.0, d5 = 0.0;

  //weight in x direction
  wx = k * (dt/(dx*dx));
  //weight in y direction
  wy = k * (dt/(dy*dy));
  
  //stability
  // wx + wy must be less than .5
  if ( wx + wy > .5)
  {
    printf("Does not reach requirements for stability\n");
    exit(1);
  }

  //MPI
  MPI_Datatype MPI_CLM;
  MPI_Type_vector(rows, 1, cols, MPI_DOUBLE, &MPI_CLM);
  MPI_Type_commit(&MPI_CLM);
  MPI_Status stat[8]; 
  MPI_Request req[8];
  int u, d, l, r;
  int startrow, startcol, endrow, endcol;

  //sets when receiver/sender process doesnt exist
  u = rank - P;
  d = rank + P;
  l = rank - 1;
  r = rank + 1;

  if (u < 0) 
    u = MPI_PROC_NULL;

  if (d >= size)
    d = MPI_PROC_NULL;

  if (rank%P == 0)
    l = MPI_PROC_NULL;
  else if (rank%P == P - 1)
    r = MPI_PROC_NULL;

  startcol = 1;
  endcol = cols - 1;
  startrow = 1;
  endrow = rows - 1;
  
  if (rank%P == 0)
    startcol = 2;
  else if (rank%P == P - 1)
    endcol = cols - 2;

  if (rank/P == 0)
    startrow = 2;
  else if (rank/P == Q - 1)
    endrow = rows - 2;

  //solving edges
  if (rank/P != 0)
  {
    //not in top blocks
    i = 1;
    for(j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
  /*     if ( d1 < fabs(mat2[i][j] - mat1[i][j])) */
  /*     { */
  /*       d1 = fabs(mat1[i][j] - mat2[i][j]); */
  /* //      printf("d5= %f\n", d5); */
  /*     }  */
    }
  }

  if (rank/P != Q - 1)
  {
    //not in bottom blocks
    i = rows - 2;
    for(j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d2 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d2 = fabs(mat1[i][j] - mat2[i][j]);
  //      printf("d5= %f\n", d5);
      }
    }
  }

  if (rank%P != 0)
  {
    //not in left blocks
    j = 1;
    for(i = startrow; i < endrow; i++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d3 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d3 = fabs(mat1[i][j] - mat2[i][j]);
  //      printf("d5= %f\n", d5);
      }
    }
  }

  if (rank%P != P-1)
  {
    //not in right blocks
    j = cols - 2;
    for(i = startrow; i < endrow; i++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d4 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d4 = fabs(mat1[i][j] - mat2[i][j]);
  //      printf("d5= %f\n", d5);
      }
    }
  }

  //halo exchange
  // Send rows up and down
  MPI_Irecv(mat2[0], cols, MPI_DOUBLE, u, 0, MPI_COMM_WORLD, &req[0]);
  MPI_Irecv(mat2[rows-1], cols, MPI_DOUBLE, d, 0, MPI_COMM_WORLD, &req[1]);
  MPI_Isend(mat2[rows-2], cols, MPI_DOUBLE, d, 0, MPI_COMM_WORLD, &req[2]);
  MPI_Isend(mat2[1], cols, MPI_DOUBLE, u, 0, MPI_COMM_WORLD, &req[3]);
  // Send columns left and right
  MPI_Irecv(&mat2[0][0], 1, MPI_CLM, l, 0, MPI_COMM_WORLD, &req[4]);
  MPI_Irecv(&mat2[0][cols-1], 1, MPI_CLM, r, 0, MPI_COMM_WORLD, &req[5]);
  MPI_Isend(&mat2[0][1], 1, MPI_CLM, l, 0, MPI_COMM_WORLD, &req[6]);
  MPI_Isend(&mat2[0][cols-2], 1, MPI_CLM, r, 0, MPI_COMM_WORLD, &req[7]);


  // perform interior calculations while send/recvs are executing
  for(i = startrow; i < endrow; i++)
  {
    for(j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
    }
  }

  MPI_Waitall(8, req, stat);


  for(i = startrow; i < endrow; i++)
    {
      for(j = startcol; j < endcol; j++)
	{
	  if ( d5 < fabs(mat2[i][j] - mat1[i][j]))
	    {
	      d5 = fabs(mat1[i][j] - mat2[i][j]);
	      printf("i: %d, j: %d, mat1: %lf, mat2: %lf\n", i, j,  mat1[i][j], mat2[i][j]);
	    }
	}
    }


  //*converge = (d5);  

  
  MPI_Allgather(&d5,1,MPI_DOUBLE,all_local_converge,1,MPI_DOUBLE,MPI_COMM_WORLD);

  *converge = all_local_converge[0];

  for(i=0;i<size;i++)
  {
    if (all_local_converge[i] > *converge)
      *converge = all_local_converge[i];
  }

  //if(!rank)
  //printf("%lf\n", d1 + d2 + d3 + d4);
  /******************
need to find max value of d5 across all of the processes, this is your "converge" value
do this by making an array all_local_converge[num_processes] and an MPI_Allgather to get all of the d5 values in the array
then find the max value in the array, set converge to be equal to this value
then all of the processes will exit at the same time and there shouldn't be a deadlock anymore

MPI_Allgather(
    &d5
    1,
    MPI_DOUBLE
    all_local_converge,
    1,
    MPI_DOUBLE
MPI_COMM_WORLD);

 ********************/

//if (rank == 0){
  for(i = 0; i < rows ; i++)
  {
    for(j = 0; j < cols ; j++)
    {
      mat1[i][j] = mat2[i][j];
    }
  }

  //  PrintGrid(mat1,nrows,ncols);

  /* double **tmp; */

  /* //swap mat1 = mat2 */
  /* tmp = mat1; */
  /* mat1 = mat2; */
  /* mat2 = tmp; */

}


int main(int argc, char *argv[])
{
  int i, j;
  double **mat1, **mat2;
  //size of grid
  int nrows = 12, ncols = 12;
  double dx = 0, dy= 0, dt, time;
  double converge = 0, epsilon;
  int iter, max_iter = 1;
  int u, d, l, r;
  clock_t start, end;

  P = 2;
  Q = 2;

  rows = 2 + nrows/Q;
  cols = 2 + ncols/P;

  //MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double *all_local_converge;
  all_local_converge = malloc(size*sizeof(double));

  // width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  epsilon = .001;
  dt = .001;
  //dt = 5;
  mat1 = calloc(rows, sizeof(double *));
  for (i = 0; i < rows; i++)
    mat1[i] = calloc(cols, sizeof(double));
	
  mat2 = calloc(rows, sizeof(double *));
  for (i = 0; i < rows; i++)
    mat2[i] = calloc(cols, sizeof(double));

  InitGrid(mat1);
  InitGrid(mat2);
  time = 0;

  PrintGrid(mat1, nrows, ncols);

  start = clock();
  for(iter = 0; iter<max_iter; iter++)
  {
    time += dt;
    Compute_TimeStep(mat1, mat2, nrows, ncols, &converge, dt, dx, dy, all_local_converge);
    
    if (converge < epsilon)
      break;
    
    /* if(!rank) */
    /*   printf("%lf\n", converge); */
  }

  //  printf("Process %d finished!\n", rank);

  MPI_Barrier(MPI_COMM_WORLD);

  end = clock();
  /* if (rank == 0) */
  /* { */
  /*   printf("num of iterations: %d, with change of %lf\n", iter, converge); */
  /*   printf("Total time: %f seconds\n", (((double) end) - start)/CLOCKS_PER_SEC); */
  /* } */

  MPI_Finalize();
}

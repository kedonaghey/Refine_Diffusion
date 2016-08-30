#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>
#include <getopt.h>

int rank, size;
int P, Q;
int rows, cols;

//initializes grid boundaries
void initGrid(double** mat)
{
  int i;
  if(rank%P == 0)
  {
    for(i= 1;i<rows-1;i++)
	mat[i][1] = 20;
  } 
  else if(rank%P == P-1)
  {
    for(i=1;i<rows-1;i++)
	mat[i][cols-2] = 30;
  }

  if(rank/P == 0)
  {
    for(i=1 ;i<cols-1;i++)
	mat[1][i] = 50;
  }
  else if(rank/P == Q-1)
  {
    for(i=1;i<cols-1;i++)
	mat[rows - 2][i] = 40;
  }
}

void printGrid(double** mat, int nrows, int ncols)
{  
  int i, j, k;

  for(k=0; k<size; k++) {
    if(rank==k) {
      printf("\nRank %d:\n", rank);
      for (i = 0; i < rows; i++) {
	for (j = 0; j < cols; j++) {
	  printf("%f ", mat[i][j]);
	}
	printf("\n");
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

   MPI_Barrier(MPI_COMM_WORLD); 
}

void computeTimestep(double*** mat1_ptr, double*** mat2_ptr, int nrows, int ncols, double* converge, double dt, double dx, double dy, double* all_local_converge)
{
  int i, j;
  double wx = 0, wy = 0;
  //diffusitivity coefficient 
  double k = 1;
  double d5 = 0.0;

  double **mat1 = *mat1_ptr;
  double **mat2 = *mat2_ptr;

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
  MPI_Type_vector(rows-2, 1, cols, MPI_DOUBLE, &MPI_CLM); 
  MPI_Type_commit(&MPI_CLM);
  MPI_Status stat[8]; 
  MPI_Request req[8];
  int u, d, l, r;
  int startrow, startcol, endrow, endcol;

  //used when rcvr/sndr process doesnt exist
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

  //solves edge points
  if (rank/P != 0)
    {
    //not in top blocks
    i = 1;
    for(j = startcol; j < endcol; j++)
      {
	mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      }
  }

  if (rank/P != Q - 1)
  {
    //not in bottom blocks
    i = rows - 2;
    for(j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
    }
  }

  if (rank%P != 0)
  {
    //not in left blocks
    j = 1;
    for(i = startrow; i < endrow; i++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
    }
  }

  if (rank%P != P-1)
  {
    //not in right blocks
    j = cols - 2;
    for(i = startrow; i < endrow; i++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
    }
  }

  //halo exchange
  // Send rows up and down
  MPI_Irecv(mat2[0],      cols, MPI_DOUBLE, u, 0, MPI_COMM_WORLD, &req[0]);
  MPI_Irecv(mat2[rows-1], cols, MPI_DOUBLE, d, 1, MPI_COMM_WORLD, &req[1]);
  MPI_Isend(mat2[rows-2], cols, MPI_DOUBLE, d, 0, MPI_COMM_WORLD, &req[2]);
  MPI_Isend(mat2[1],      cols, MPI_DOUBLE, u, 1, MPI_COMM_WORLD, &req[3]);
  // Send columns left and right
  MPI_Irecv(&mat2[1][cols-1], 1, MPI_CLM, r, 2, MPI_COMM_WORLD, &req[4]); 
  MPI_Irecv(&mat2[1][0],      1, MPI_CLM, l, 3, MPI_COMM_WORLD, &req[5]); 
  MPI_Isend(&mat2[1][1],      1, MPI_CLM, l, 2, MPI_COMM_WORLD, &req[6]); 
  MPI_Isend(&mat2[1][cols-2], 1, MPI_CLM, r, 3, MPI_COMM_WORLD, &req[7]); 


  // perform interior calculations while send/recvs are running
  for(i = startrow; i < endrow; i++)
  {
    for(j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
    }
  }

  MPI_Waitall(8, req, stat);

  //calculates the converge value
  for(i = startrow; i < endrow; i++)
    {
      for(j = startcol; j < endcol; j++)
	{
	  if ( d5 < fabs(mat2[i][j] - mat1[i][j]))
	    {
	      d5 = fabs(mat1[i][j] - mat2[i][j]);
	    }
	}
    }

  // all processes need to know the convergence value
  MPI_Allgather(&d5, 1, MPI_DOUBLE, all_local_converge, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  *converge = all_local_converge[0];
  //finds the correct value
  for(i=0;i<size;i++)
  {
    if (all_local_converge[i] > *converge)
      *converge = all_local_converge[i];
  }


  //swap mat1 = mat2
  double** tmp = *mat1_ptr;
  *mat1_ptr = *mat2_ptr;
  *mat2_ptr = tmp; 
}


int main(int argc, char *argv[])
{
  int i, c = 0;
  double **mat1, **mat2;
  //size of grid
  int nrows = 12, ncols = 12;
  double dx = 0, dy= 0, dt, time;
  double converge = 0, epsilon;
  int iter, max_iter = 1000;
  clock_t start, end;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  //P * Q = num procs
  P = 2;
  Q = 2;


  while ((c = getopt (argc, argv, "p:q:r:i:")) != -1)
  {
    switch(c)
    {
      case 'p':
	P = ncols = atoi(optarg); break;
      case 'q':
	Q = atoi(optarg); break;
      case 'r':
	nrows = ncols = atoi(optarg); break;
      case 'i':
	max_iter = atoi(optarg); break;
      default:
 	fprintf(stderr, "Invalid option\n");
	return -1;
    }
  }

  //if size doesnt equal P * Q print error message
  if(size != P*Q)
  {
    if (!rank)
	printf("P*Q must run with %d MPI tasks\n", P*Q);
    MPI_Finalize();
    exit(0);
  }

  //if not divisble by P or Q, print error message
  if(nrows%Q != 0 || nrows%P != 0)
  {
    if (!rank)
      printf("P %d and Q %d must divide evenly into %d \n", P, Q, nrows);
    MPI_Finalize();
    exit(0);
  }
  
  //adds halo rows
  rows = 2 + nrows/Q;
  cols = 2 + ncols/P;

  double *all_local_converge;
  all_local_converge = malloc(size*sizeof(double));

  // width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  epsilon = .001;
  dt = .00001;

  // allocate matrices in 1d array with 2d access
  mat1 = calloc(rows, sizeof(double *)); 
  mat2 = calloc(rows, sizeof(double *));
  double *mat1_1d = calloc(rows*cols, sizeof(double));
  for(i=0; i<rows; i++) {
    mat1[i] = &(mat1_1d[i*cols]);
  }

  double *mat2_1d = calloc(rows*cols, sizeof(double));
  for(i=0; i<rows; i++) {
    mat2[i] = &(mat2_1d[i*cols]);
  }

  initGrid(mat1);
  initGrid(mat2);
  time = 0;
 //iterates through until convergence is met
 start = clock();
 for(iter = 0; iter<max_iter; iter++)
 {
   time += dt;
   computeTimestep(&mat1, &mat2, nrows, ncols, &converge, dt, dx, dy, all_local_converge);
   
   if (converge < epsilon)
     break;
   
     if(!rank)  
       printf("iter %d: %lf\n",iter, converge);  
 }

 MPI_Barrier(MPI_COMM_WORLD);

 end = clock();
   if (rank == 0)  
   {  
     printf("num of iterations: %d, with change of %lf\n", iter, converge);  
     printf("Total time: %f seconds\n", (((double) end) - start)/CLOCKS_PER_SEC);  
   }  

  MPI_Finalize();

  return 0;
}

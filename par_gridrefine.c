#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

int rank, size, cRank, cRank2;
MPI_Comm MPI_COMM_CART;
MPI_Comm MPI_COMM_CART2;
int P, Q, R, S;
int crds[2];
int rcrds[2];
int par_rows, par_clms;
int par_rrows, par_rclms;
MPI_Datatype MPI_CLM;
MPI_Datatype MPI_ROW;

void initFineGrid(double** mat, int nrows, int ncols)
{
  int i;
  int top = 20;
  int left = 50;
  //bottom and top rows initial conditions
  if(rcrds[1] == R-1)
  {
  for ( i = 0; i < nrows; i+=2)
    {
      mat[i][0] = top;
    }
  for (i = 1; i < nrows; i+=2)
    {
      mat[i][0] = top/2;
    }
  }
  //left and right rows initial conditions
  if(rcrds[0] == 0)
  {
  for ( i = 0; i < ncols; i+=2)
    {
      mat[nrows-1][i] = left;//100
    }
  for (i=1; i < ncols; i+=2)
    {
      mat[nrows-1][i] = left/2;
    }
  } 
}

void initGrid(double** mat, int nrows, int ncols, int crs, int ccs)
{
  int i;
  //bottom and top rows initial conditions
  if(crds[1] == 0)
  {
  for ( i = 0; i < nrows - crs +1; i++)
    {
      mat[i][0] = 20;//0
    }
  }
  else if(crds[1] == P - 1)
  {
  for (i = 0; i <nrows; i++)
    { 
      mat[i][ncols - 1] = 30;//100
    }
  }

  //left and right rows initial conditions
  if(crds[0] == 0)
  {
  for ( i = ccs -1; i < ncols; i++)
    {
      mat[nrows - 1][i] = 50;//0
    }
  }
  else if(crds[0] == Q-1)
  {
  for ( i = 0; i < ncols; i++)
    {
      mat[0][i] = 40;//100
    }
  }

}

void printGrid(double** mat, int nrows, int ncols)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }


}
void printGridToFile(double** mat, int nrows, int ncols, FILE *fp)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            fprintf(fp, "%f ", mat[i][j]);
        }
        fprintf(fp, "\n");
    }


}

void injectCoarsePoints(double* toparray, double* rightarray, double** mat, int nrows, int ncols, int crs, int ccs, double corner)
{
  int i, j = 0;
  //int test[] = {1,2,3};
  //top
  for (i = 1; i < ccs; i++)
      mat[nrows-crs][i] = toparray[i-1];//toparray[i-1];//0
  //right
  for(i = nrows-2; i > nrows - crs; i--)
  {
      //mat[i][ccs - 1] = test[i-crs-1];//rightarray[i-crs-1];//100
      mat[i][ccs - 1] = rightarray[j];
      j++;
  }
  //corner
  mat[nrows-crs][ccs-1] = corner;
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

  mat2_corner = /*prev pointmat1_corner*/mat1[i][j] + wxy * (16/15) * (/*prev point*/-4* mat1[i][j]/*mat1_corner*/+\
  .5 * /*bound lwr*/mat1[i+1][j] + /*coarse*/mat1[i][j+1] + /*bound*/ mat1[i-1][j] + .5 * /*coarse*/mat1[i][j-1]\
  +/*mesh*/ mat1_refine[m+1][n-1]);

 return mat2_corner;
}

//need to add in different

void computeInterfaceRightTimestep(double* rightarray, double** mat1_refine, double** mat1, int nrows, int ncols, int ccs, int m, double wxy, double dx, int crs, int refcols)
{
  int j = ccs-1, n = refcols-1, t = 0, i;
  m = 2;

    for (i = nrows-2; i > nrows-crs; i--)
    {
      rightarray[t] =/*previous interface point*/ mat1[i][j] + /*previous interface point*/ mat1[i][j] * wxy * (dx/2) + wxy *\
      (4/3) * (/*previous interface point*/- 4*mat1[i][j] + /*coarse*/mat1[i][j+1]  + /*coarse*/ .5 * mat1[i+1][j]\
      +/*coarse*/ .5 * mat1[i-1][j] +/*refine*/ .5 * mat1_refine[m-1][n-1] +/*refine*/ mat1_refine[m][n-1] +/*refine*/ .5\
      * mat1_refine[m+1][n-1]);//changed one n to - different from the book
      m+=2;
      //iterates through right array
      t++;
    }


}


void computeInterfaceTopTimestep(double* toparray, double** mat1_refine, double** mat1, int nrows, int ncols, int crs, int n, double wxy, double dx, int ccs)
{
  int i = nrows - crs;
  int m = 0, t = 0, j;
  n = 2;

  for (j = 1; j < ccs; j++)
  {
      toparray[t] = /*previous interface point*/mat1[i][j] + /*prev interface point*/mat1[i][j] * wxy * (dx/2) + wxy * (4/3) * \
      (-/*prev interface point*/ 4*mat1_refine[i][j] +/*coarse*/ mat1[i-1][j]  +/*coarse*/ .5 * mat1_refine[i][j+1]\
      +/*coarse*/ .5 * mat1_refine[i][j-1] +/*refine*/ .5 * mat1_refine[m+1][n-1] +/*refine*/ mat1_refine[m+1][n] +/*refine*/ .5 * mat1_refine[m+1][n+1]);
      n+=2;
      //iterates through top array
      t++;
  }

}


void computeFineTimestep(double** mat1_refine, double** mat2_refine, int nrows, int ncols, double* converge, double wxy, double sxy)
{

  int i,j;
  MPI_Status stat2[8]; 
  MPI_Request req2[8];
  int u, d, l, r;
  int startrow, startcol, endrow, endcol;

  MPI_Cart_shift(MPI_COMM_CART, 0, 1, &u, &d);
  MPI_Cart_shift(MPI_COMM_CART, 1, 1, &l, &r);

  startcol = 1;
  endcol = par_rcols - 1;
  startrow = 1;
  endrow = par_rrows - 1;

  if (rcrds[1] == 0)
    startcol = 2;
  else if (rcrds[1] == P - 1)
    endcol = par_rcols -2;

  if (rcrds[0] == 0)
    startrow = 2;
  else if (rcrds[0] == Q - 1)
    endrow = par_rrows -2;

  //solve edges
  if (rcrds[0] != 0)
  {
  //not in top blocks
  i = 1;
    for (j = startcol; j < endcol; j++)
    {
      if (j%2 == 1)
      {
      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i-1][j+1]\
      + mat1_refine[i-1][j-1]) -2 * mat1_refine[i][j] + mat1_refine[i+1][j]) + wxy * (mat1_refine[i][j-1] -2 * mat1_refine[i][j]\
      + mat1_refine[i][j+1]);
      }
      else
      {
      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]\
      -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
      }
    }
  }
      
  if (rcrds[0] != Q-1)
  {
  //not in bottom blocks
  i = par_rows - 2;
  for (j = 1; j < ncols - 1; j++)
    {
      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]\
      -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
    }
  }

  if (rcrds[1] != 0)
  {
  j = 1;
  for (i = startrow; i < endrow; i++)
      {
      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * (mat1_refine[i+1][j]\
      -2 * mat1_refine[i][j] + mat1_refine[i-1][j]) + wxy * (mat1_refine[i][j+1] - 2 * mat1_refine[i][j] + mat1_refine[i][j-1]);
      }
  }


  if (rcrds[1] != 0)
  {
  j = par_cols - 2;
  for (i = startrow; i < endrow; i++)
  {
      if (i%2 == 1)
      {
      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i+1][j+1]\
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

  //halo exchange
  // Send rows up and down
  MPI_Irecv(mat2[0],      par_rcols, MPI_DOUBLE, u, 0, MPI_COMM_CART2, &req2[0]);
  MPI_Irecv(mat2[par_rows-1], par_rcols, MPI_DOUBLE, d, 1, MPI_COMM_CART2, &req2[1]);
  MPI_Isend(mat2[par_rows-2], par_rcols, MPI_DOUBLE, d, 0, MPI_COMM_CART2, &req2[2]);
  MPI_Isend(mat2[1],      par_rcols, MPI_DOUBLE, u, 1, MPI_COMM_CART2, &req2[3]);
  // Send columns left and right
  MPI_Irecv(&mat2[1][par_rcols-1], 1, MPI_CLM, r, 2, MPI_COMM_CART2, &req2[4]);
  MPI_Irecv(&mat2[1][0],      1, MPI_CLM, l, 3, MPI_COMM_CART2, &req2[5]);
  MPI_Isend(&mat2[1][1],      1, MPI_CLM, l, 2, MPI_COMM_CART2, &req2[6]);
  MPI_Isend(&mat2[1][par_rcols-2], 1, MPI_CLM, r, 3, MPI_COMM_CART2, &req2[7]);

  for (i = startrow; i < endrow; i++)
  {
    for (j = startrow; j < endrow; j++)
    {
      if (i == 1 && j%2 == 1)
      {
      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i-1][j+1]\
      + mat1_refine[i-1][j-1]) -2 * mat1_refine[i][j] + mat1_refine[i+1][j]) + wxy * (mat1_refine[i][j-1] -2 * mat1_refine[i][j]\
      + mat1_refine[i][j+1]);
      }
      if (j == ncols-2 && i%2 == 1)
      {
      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i+1][j+1]\
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


void computeTimestep(double** mat1, double** mat2, int nrows, int ncols, double* converge, double wx, double wy, int crs, int ccs)
{
  int i, j;
  MPI_Status stat[8]; 
  MPI_Request req[8];
  int u, d, l, r;
  int startrow, startcol, endrow, endcol;

  MPI_Cart_shift(MPI_COMM_CART, 0, 1, &u, &d);
  MPI_Cart_shift(MPI_COMM_CART, 1, 1, &l, &r);

  startcol = 1;
  endcol = par_cols - 1;
  startrow = 1;
  endrow = par_rows - 1;

  if (crds[1] == 0)
    startcol = 2;
  else if (crds[1] == P - 1)
    endcol = par_cols -2;

  if (crds[0] == 0)
    startrow = 2;
  else if (crds[0] == Q - 1)
    endrow = par_rows -2;

  //solve edges
  if (crds[0] != 0)
  {
  //not in top blocks
  i = 1;
    for (j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
}}

  if (crds[0] != Q-1)
  {
  //not in bottom blocks
  i = par_rows - 2;
  for (j = 1; j < ncols - 1; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
}}

  if (crds[1] != 0)
  {
  j = 1;
  for (i = startrow; i < endrow; i++)
  {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
}}

  if (crds[1] != 0)
  {
  j = par_cols - 2;
  for (i = startrow; i < endrow; i++)
  {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
}}

  //halo exchange
  // Send rows up and down
  MPI_Irecv(mat2[0],      cols, MPI_DOUBLE, u, 0, MPI_COMM_CART, &req[0]);
  MPI_Irecv(mat2[par_rows-1], cols, MPI_DOUBLE, d, 1, MPI_COMM_CART, &req[1]);
  MPI_Isend(mat2[par_rows-2], cols, MPI_DOUBLE, d, 0, MPI_COMM_CART, &req[2]);
  MPI_Isend(mat2[1],      cols, MPI_DOUBLE, u, 1, MPI_COMM_CART, &req[3]);
  // Send columns left and right
  MPI_Irecv(&mat2[1][par_cols-1], 1, MPI_CLM, r, 2, MPI_COMM_CART, &req[4]);
  MPI_Irecv(&mat2[1][0],      1, MPI_CLM, l, 3, MPI_COMM_CART, &req[5]);
  MPI_Isend(&mat2[1][1],      1, MPI_CLM, l, 2, MPI_COMM_CART, &req[6]);
  MPI_Isend(&mat2[1][par_cols-2], 1, MPI_CLM, r, 3, MPI_COMM_CART, &req[7]);


  for (i = startrow; i < endrow; i++)
  {
    for (j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
}}

  MPI_Waitall(8,req, stat);

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

  //MPI Init
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  P = 2;
  Q = 2;

  R = 2;
  S = 2;


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

  //parallel rows & clms coarse grid
  int par_rows = 2 + nrows/Q;
  int par_clms = 2 + ncols/P;

  int szofmesh = .4*nrows;

  //size of fine grid
  int refrows = szofmesh * 2 + 1;
  int refcols = szofmesh * 2 + 1;

  //parallel refined rows & clms
  int par_rrows = 2 + refrows/R;
  int par_rcols = 2 + refcols/S;


  MPI_Type_vector(par_rrows - 2, 1, par_rcols * 2, MPI_DOUBLE, &MPI_CLM);
  MPI_Type_commit(&MPI_CLM);

//  MPI_Type_vector(par_crows - 2, 1, par_cols * 2, MPI_DOUBLE, &MPI_CLM);
//  MPI_Type_commit(&MPI_ROW);

  //Create cart
  int prds[] = {0,0};
  int dims[] = {Q,P};
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, prds, 1, &MPI_COMM_CART);
  MPI_Comm_rank(MPI_COMM_CART, &cRank);
  MPI_Cart_coords(MPI_COMM_CART, cRank, 2, crds);

  //Create cart
  int rprds[] = {0,0};
  int rdims[] = {R,S};
  MPI_Cart_create(MPI_COMM_WORLD, 2, rdims, rprds, 1, &MPI_COMM_CART2);
  MPI_Comm_rank(MPI_COMM_CART2, &cRank2);
  MPI_Cart_coords(MPI_COMM_CART2, cRank2, 2, rcrds);


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

  int crs = refrows/2 + 1;
  int ccs = refcols/2 + 1;
  initGrid(mat1, nrows, ncols, crs, ccs);
  initGrid(mat2, nrows, ncols, crs, ccs);

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
  int ci = nrows - crs;
  int cj = ccs - 1;
  //fine
  int fi = 0;
  int fj = refcols - 1;

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
      computeTimestep(mat1, mat2, nrows, ncols, &converge, wx, wy, crs, ccs);
      computeFineTimestep(mat1_refine, mat2_refine, refrows, refcols, &converge, wxy, sxy);
      //Compute Interfaces
      computeInterfaceRightTimestep(rightarray, mat1_refine, mat1, nrows, ncols, ccs, fi, wxy, dx, crs, refcols);
      computeInterfaceTopTimestep(toparray, mat1_refine, mat1, nrows, ncols, crs, fj, wxy, dx, ccs);
      double corner = computeCornerTimestep(mat1_refine, mat1, ci, cj, fi, fj, wxy);
      injectCoarsePoints(toparray, rightarray, mat2, nrows, ncols, crs, ccs, corner);
      injectFinePoints(toparray, rightarray, mat2_refine, refrows, refcols, corner);
      update(&mat1, &mat2, nrows, ncols);
      update(&mat1_refine, &mat2_refine, refrows, refcols);
      //if (iter%100 == 0)
      iter++;
    }
  //end = clock();
  //for(i = 0; i < szofmesh -1; i++)
  //printf("mat1 %lf\n", mat1[nrows-crs][3-1]);
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


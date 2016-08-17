#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#define PI 3.14159265358979
//top is 100 everything else 0
int main(int argc,char* argv[])
{
  int m, n, i, j;
  double** mat;
  double** finmat;
  double lambda, B, k = 1; //diffusion coefficient
  int t = atoi(argv[1]);
  double mn;
  //point on grid in x direction
  int x,
  //point on grid in y direction
  y,
  //total length in y direction
  y0 = 10,
  //total length in x direction
  x0 = 10;

  double c;

  //max num of iterations
  int max_value = 100;

  mat = calloc(x0, sizeof(double *));
  for (i = 0; i < x0; i++)
    mat[i] = calloc(y0, sizeof(double));
  
  for (x = 0; x < x0; x++)
  {
  for (y = 0; y < y0; y++)
  {
  for(m = 1; m < max_value; m++)
  {
    for(n = 1; n < max_value; n++)
    {
      c = (200/(PI*PI)) * (((cos(PI*m)-1)*(cos(PI*n)-1))/(m*n));
      lambda = ((double)(m*m)/(x0*x0) + (double)(n*n)/(y0*y0)) * (PI*PI*k);
      mat[x][y] +=  c * sin((PI*x*m)/x0) * sin((PI*y*n)/(y0)) * exp(-1*lambda*t);
    }
  }
     
  }}

//when m and n = 1

/*
  for (x = 0; x < x0; x++)
  {
    for (y = 0; y < y0; y++)
    {
      c = (200/(PI*PI)) * (((cos(PI*1)-1)*(cos(PI*1)-1))/(1*1));
      lambda = ((1*1)/(x0*x0) + (1*1)/(y0*y0)) * (PI*PI*k);
      mat[x][y] +=  c * sin((PI*x*1)/x0) * sin((PI*y*1)/(y0)) * exp(-1*lambda*t);
    }
  }
*/
  for (i = 0; i < x0; i++) 
  {
    for (j = 0; j < y0; j++) 
    {
      printf("%f ", mat[i][j]);
    }
    printf("\n");
  }


}




/*

      //printf("m = %d, n = %d, x0 = %d, y0 =%d, PI = %lf, k = %lf, lambda = %20.16lf\n",m,n,x0,y0,PI,k,lambda);
      //printf("m = %lf\n", (double)(m*m)/(x0*x0);
      mn = c * sin((PI*x*m)/x0) * sin((PI*y*n)/(y0));
      //if (x==5 && y==5)
        //printf("m = %d, n =%d; lambda = %lf, mn = %lf\n", m,n, lambda, mn);
*/

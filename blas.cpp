/*---------------------------------------------------------------------------+
|  g++ -o blas blas.cpp -lblas                                               |
+---------------------------------------------------------------------------*/
#define for_all_x(i) for(int i=0; i<n_x; i++)
#define for_all_y(j) for(int j=0; j<n_y; j++)
#define for_all_z(k) for(int k=0; k<n_z; k++)
#define for_all_xyz(i,j,k) for_all_x(i) for_all_y(j) for_all_z(k)

extern "C" {
  extern void daxpy_(int * n, double * a, double * x, int * incx, double * y, int * incy);
};
void allocate(double **** val, const int n_x, const int n_y, const int n_z);

void daxpy(int n, double a, double * x, double * y) {
  for(int i=0; i<n; i++)
    y[i] += a * x[i];
}


/******************************************************************************/
#include <iostream>
int main(int argc, char * argv[]) {

  if(argc != 2) {
    std::cout << "Right syntax: a.out [3 / L / B]" << std::endl;
    return -1;
  }

  const int n_x=100, n_y=100, n_z=200;

  double a = 10;
  double *** x, *** y;

  allocate(&x, n_x, n_y, n_z);
  allocate(&y, n_x, n_y, n_z);

  for_all_xyz(i,j,k) {
    x[i][j][k] = 10;
    y[i][j][k] =  1;
  }

  /*--------------------------+
  |  Threedimensional access  |
  +--------------------------*/
  if( argv[1][0] == '3' ) {
    std::cout << "Testing 3D access" << std::endl;
    for(int k=0; k<500; k++)
      for_all_xyz(i,j,k)
        y[i][j][k] += a * x[i][j][k];
  }

  /*----------------+
  |  Linear access  |
  +----------------*/
  if( argv[1][0] == 'L' ) {
    std::cout << "Testing linear access" << std::endl;
    int n_t = n_x*n_y*n_z;
    for(int k=0; k<500; k++)
      daxpy(n_t, a, x[0][0], y[0][0]);
  }

  /*-------+
  |  BLAS  |
  +-------*/
  if( argv[1][0] == 'B' ) {
    std::cout << "Testing BLAS" << std::endl;
    int inc = 1;
    int n = n_x * n_y * n_z;
    for(int k=0; k<500; k++)
      daxpy_(&n, &a, x[0][0], &inc, y[0][0], &inc);
  }

  std::cout << y[0][0][0] << std::endl;
  std::cout << y[50][50][50] << std::endl;
}

/******************************************************************************/
void allocate(double **** val, const int n_x, const int n_y, const int n_z) {

  /* allocate memory as a contiguous block */
  (*val)    = new double ** [n_x];

  (*val)[0] = new double *  [n_x * n_y];
  for_all_x(i)
    (*val)[i] = (*val)[0] + i*n_y;

  (*val)[0][0] = new double [n_x * n_y * n_z]; // important: contiguous block
  for_all_x(i)
    for_all_y(j)
      (*val)[i][j] = (*val)[0][0] + i*n_y*n_z + j*n_z;

  std::cout << "Allocation OK" << std::endl;
}

/*-----------------------------------------------------------------------------+
 '$Id: blas.cpp,v 1.4 2008/11/17 19:23:21 niceno Exp $'/
+-----------------------------------------------------------------------------*/

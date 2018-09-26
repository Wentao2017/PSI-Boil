/*----------------------------------------------------------------------------+
|  CC -o blas_2 -O4 -fastsse -Minfo=loop -Mneginfo=loop -Msafeptr blas_2.cpp  |
+----------------------------------------------------------------------------*/

extern "C" {
  extern void saxpy_(int * n, float * a, float * x, int * incx, float * y, int * incy);
};
void allocate(float **** val, const int NX, const int NY, const int NZ);

/******************************************************************************/
#include <iostream>
int main(int argc, char * argv[]) {

  if(argc != 2) {
    std::cout << "Right syntax: a.out [N / B]" << std::endl;
    return -1;
  }

  const int NX=256, NY=256, NZ=256;

  float a = 10;
  float *** x, *** y;

  allocate(&x, NX, NY, NZ);
  allocate(&y, NX, NY, NZ);

  /*-------------+
  |  Initialize  |
  +-------------*/
  for(int i=0; i<NX; i++)  
    for(int j=0; j<NY; j++)
      for(int k=0; k<NZ; k++) { 
    x[i][j][k] = 10;
    y[i][j][k] =  1;
  }

  for(int k=0; k<512; k++) {

    /*----------------+
    |  Normal access  |
    +----------------*/
    if( argv[1][0] == 'N' ) {
      std::cout << "Testing normal access" << std::endl;
        for(int i=0; i<NX; i++)
          for(int j=0; j<NY; j++)
            for(int k=0; k<NZ; k++) 
          y[i][j][k] += a * x[i][j][k];
    }

    /*-------+
    |  BLAS  |
    +-------*/
    if( argv[1][0] == 'B' ) {
      std::cout << "Testing BLAS" << std::endl;
      int inc = 1;
      int n = NX * NY * NZ;
        saxpy_(&n, &a, x[0][0], &inc, y[0][0], &inc);
    }
  }

}

/******************************************************************************/
void allocate(float **** val, const int NX, const int NY, const int NZ) {

  /* allocate memory as a contiguous block */
  (*val)    = new float ** [NX];

  (*val)[0] = new float *  [NX * NY];
  for(int i=0; i<NX; i++)
    (*val)[i] = (*val)[0] + i*NY;

  (*val)[0][0] = new float [NX * NY * NZ]; // important: contiguous block
  for(int i=0; i<NX; i++)
    for(int j=0; j<NY; j++)
      (*val)[i][j] = (*val)[0][0] + i*NY*NZ + j*NZ;

  std::cout << "Allocation OK" << std::endl;
}

/*-----------------------------------------------------------------------------+
 '$Id: blas_2.cpp,v 1.5 2008/11/29 13:39:54 niceno Exp $'/
+-----------------------------------------------------------------------------*/

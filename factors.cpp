#include <iostream>

using namespace std;

/******************************************************************************/
void swap(int & a, int & b) {
  int t = b;
  b = a;
  a = t;
}

/******************************************************************************/
void factor(int n, int * factor, int * number) {
 
  (*number) = 0;
  
  int m=n-1;
  for(;;) {
    if(n % m == 0) {
      factor[(*number)++] = n/m;
      n = m;
    }
    m--;	  
    if(m == 0) break;
    if(m == 1) {
      factor[(*number)++] = n;
      break;
    }
  }
}

/******************************************************************************/
void distribute(const int nproc, const int dim, int * dis, int * res) {

  int factors[64], nfact;
  factor( nproc, factors, &nfact);      	
  
  /*------------------+
  |                   |
  |  initializations  |
  |                   |
  +------------------*/
  for(int m=0; m<dim; m++)
    dis[m] = 1;
  
  /* actual dimensions (not taking the freezed ones) */
  int act_dim = dim;
  for(int m=0; m<dim; m++)
    if(res[m] == 1) act_dim--;

  /*--------------------------------+
  |                                 |
  |  define unsorted distributions  |
  |                                 |
  +--------------------------------*/

  /*------------------------------------------------------+
  |  number of dimensions equal to the number of factors  |
  +------------------------------------------------------*/
  if(act_dim == nfact) {
    for(int m=0; m<act_dim; m++)
      dis[m] = factors[m];
  }

  /*-------------------------------+
  |  less factors than dimensions  |
  +-------------------------------*/
  if(act_dim > nfact) {
    for(int m=0; m<nfact; m++)
      dis[m] = factors[m];
  }

  /*-------------------------------+
  |  more factors than dimensions  |
  +-------------------------------*/
  if(act_dim < nfact) {
	  
    /* backward first */
    int l = nfact-1;	  
    for(int m=act_dim-1; m>=0; m--)
      dis[m] *= factors[l--];

    /* then forward */
    int m=0;
    for(int l=0; l<nfact-act_dim; l++) {
      dis[m++] *= factors[l];
      if(m==act_dim) m=0;
    }
  }

  /*------------------------------+
  |                               |
  |  (bubble) sort distributions  |
  |                               |
  +------------------------------*/
  int sor_res[] = {res[0], res[1], res[2]};
  int sor_dis[] = {dis[0], dis[1], dis[2]};
  int pos_res[] = { 0, 1, 2};
  for(int n=0; n<dim; n++) {
    for(int m=0; m<dim-1; m++) {

      if(sor_dis[m] > sor_dis[m+1]) {swap(sor_dis[m], sor_dis[m+1]);}
      if(sor_res[m] > sor_res[m+1]) {swap(sor_res[m], sor_res[m+1]);
                                     swap(pos_res[m], pos_res[m+1]);}
    }
  }

  for(int m=0; m<dim; m++) 
    dis[ pos_res[m] ] = sor_dis[ m ];
}
 
/******************************************************************************/
int main(int argc, char ** argv) {
 
  /*-------------------+
  |  check the sintax  |
  +-------------------*/
  if(argc < 4 || argc > 5 ) {
    cout << "Syntax: factor Nproc Nx Ny [Nz]" << endl;
    return -1;
  }

  int nproc  = atoi(argv[1]);
  int ndim   = argc - 2;

  /* resolution in each direction */
  int res[3];
  res[0] = atoi(argv[2]);
  res[1] = atoi(argv[3]);
  res[2] = 1;
  if(argc == 5)
    res[2] = atoi(argv[4]);
  
  /* distribution over processors */
  int dis[3];

  distribute(nproc, ndim, dis, res);
  
  for(int m=0; m<ndim; m++)
    cout << "res|dis = " << res[m] << "|" << dis[m] << endl;

  return 0;
}

/*-----------------------------------------------------------------------------+
 '$Id: factors.cpp,v 1.5 2008/11/29 13:39:54 niceno Exp $'/
+-----------------------------------------------------------------------------*/

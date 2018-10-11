#include "vof.h"
#include <cmath>
#include <fstream>

/******************************************************************************/
void VOF::plic() {
#if 0
  const real dt= time->dt();
  real dx, dy, dz;
  real dtdx, dtdy, dtdz;
  Comp mc;
#endif
 
  std::cout<<"si(x) "<<u->si(Comp::u())<<"ei(x) "<<u->ei(Comp::u())<<"si "<<si()<<"ei "<<ei()<<"\n"; 

  gradphic(phi);
  plane_vector_mc();
  calc_alpha();

#if 0
  /*-------------------+
  |  Backup variables  |
  +--------------------*/

  for_ijk(i,j,k){

    alpha_tmp[i][j][k] = alpha[i][j][k];
    vma_tmp[i][j][k]   = vma[i][j][k];
    vmb_tmp[i][j][k]   = vmb[i][j][k];
    vmc_tmp[i][j][k]   = vmc[i][j][k];

  }

  /*----------------------+
  |  Flux in X direction  |
  +----------------------*/

  for_ijk(i,j,k){

    dx   = phi.dxc(i);
    dtdx = dt/dx;

    mc = Comp::u();
    absgu[i][j][k] = fabs((*u)[mc][i][j][k] * dtdx);
   
    ra[i][j][k] = vma[i][j][k] * (1.0 - absgu[i][j][k]);
    qa[i][j][k] = 1.0 / (1.0 - ra[i][j][k]);

    if((*u)[mc][i][j][k] * dtdx > 0.0){
      
      alpha[i][j][k] = alpha[i][j][k] - ra[i][j][k];

    }
    
    vma[i][j][k]   = vma[i][j][k] * absgu[i][j][k];  

    alpha[i][j][k] = alpha[i][j][k] * qa[i][j][k];
    vma[i][j][k]   = vma[i][j][k]   * qa[i][j][k]; 
    vmb[i][j][k]   = vmb[i][j][k]   * qa[i][j][k];
    vmc[i][j][k]   = vmc[i][j][k]   * qa[i][j][k];

  }

  calc_v();

  for_ijk(i,j,k){

    flux_x[i][j][k] = vv[i][j][k] * (*u)[mc][i][j][k] * dtdx;
  
  }

  /*----------------------+
  |  Flux in Y direction  |
  +----------------------*/

  for_ijk(i,j,k){

    alpha[i][j][k] = alpha_tmp[i][j][k];
    vma[i][j][k]   = vma_tmp[i][j][k];
    vmb[i][j][k]   = vmb_tmp[i][j][k];
    vmc[i][j][k]   = vmc_tmp[i][j][k];
 
  }

  for_ijk(i,j,k){

    dy    =   phi.dyc(j);
    dtdy  =   dt/dy;

    mc = Comp::v();
    absgv[i][j][k] = fabs((*v)[mc][i][j][k] * dtdy); 

    ra[i][j][k] = vmb[i][j][k] * (1.0 - absgv[i][j][k]); 
    qa[i][j][k] = 1.0 / (1.0 - ra[i][j][k]);


#if 1
  boil::plot->plot(phi,nx,ny,nz, "clr-nx-ny-nz");
  boil::plot->plot(phi,vma,vmb,vmc, "clr-vma-vmb-vmc");
  boil::plot->plot(alpha, "alpha");
#endif

#endif
  return;
}
  

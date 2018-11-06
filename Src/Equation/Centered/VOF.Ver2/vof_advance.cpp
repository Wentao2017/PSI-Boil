#include "vof.h"

/******************************************************************************/
void VOF::advance() {

  gradphic(phi);

  // advance in x-direction
   advance_x();


  // advance in y-direction
   advance_y();
  

  // advance in z-direction
 //  advance_z();

 
  flx.exchange();
  
  for_ijk(i,j,k){
    real phi_tmp = phi[i][j][k] * dV(i,j,k);
    phi_tmp +=  flx[Comp::u()][i][j][k] - flx[Comp::u()][i+1][j][k]
              + flx[Comp::v()][i][j][k] - flx[Comp::v()][i][j+1][k]
              + flx[Comp::w()][i][j][k] - flx[Comp::w()][i][j][k+1];
    phi_tmp /= dV(i,j,k);
    phi[i][j][k] = std::min(1.0,std::max(0.0,phi_tmp));

#if 1
    if(phi_tmp>1.0+boil::pico || phi_tmp< -boil::pico){
      std::cout.setf(std::ios_base::scientific);
      std::cout<<"limit phi "<<phi_tmp<<" i "<<i<<" j "<<j<<" k "<<k<<"\n";
      std::cout.unsetf(std::ios_base::floatfield);
    }
#endif

  }
  phi.bnd_update();
  phi.exchange_all();

  return;
}


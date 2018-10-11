#include "vof.h"

/******************************************************************************/
void VOF::advance() {

  stmp = phi;

  // advance in x-direction
   advance_x();


  // advance in y-direction
   advance_y();
  

  // advance in z-direction
 //  advance_z();

  // limit phi between 0 <= phi <= 1.0
  for_ijk(i,j,k){
    real phi_tmp = phi[i][j][k];
#if 1
    if(phi_tmp>1.0+boil::pico || phi_tmp< -boil::pico){
      std::cout.setf(std::ios_base::scientific);
      std::cout<<"limit phi "<<phi_tmp<<" i "<<i<<" j "<<j<<" k "<<k<<"\n";
      std::cout.unsetf(std::ios_base::floatfield);
    }
#endif
    phi[i][j][k] = std::min(1.0,std::max(0.0,phi_tmp));
  }
  phi.bnd_update();
  phi.exchange_all();

  return;
}


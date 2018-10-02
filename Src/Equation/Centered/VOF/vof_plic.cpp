#include "vof.h"
#include <cmath>

/******************************************************************************/
void VOF::plic() {
  
  gradphic(phi);
  plane_vector_mc();

#if 1
  boil::plot->plot(phi,nx,ny,nz, "clr-nx-ny-nz");
  boil::plot->plot(phi,vma,vmb,vmc, "clr-vma-vmb-vmc");
#endif

  return;
}
  

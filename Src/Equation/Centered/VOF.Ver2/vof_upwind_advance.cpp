#include "vof.h"
#include <cmath>

/******************************************************************************/
void VOF::upwind_advance() {

phi.bnd_update();
phi.exchange_all();

const real dt= time->dt();
real dx, dy, dz;
real dtdx, dtdy, dtdz;
real dx2, dx3, dy2, dy3, dz2, dz3;
real cx, cy, cz;
real fluxx, fluxy,fluxz;
Comp mc;

    std::cout<< phi[80][1][1] <<"\n";
#if 1
for_ijk(i,j,k){

    dx=phi.dxc(i);
    dy=phi.dyc(j);
    dz=phi.dzc(k);

    dtdx=dt/dx;
    dtdy=dt/dy;
    dtdz=dt/dz;

    dx2 =dx*dx;
    dx3 =dx2*dx;
    dy2 =dy*dy;
    dy3 =dy2*dy;
    dz2 =dz*dz;
    dz3 =dz2*dz;

  // std::cout<< clr[80][1][1] <<"\n";
/* Flux */
    mc = Comp::u();
    fluxx = phi[i][j][k] * (*u)[mc][i+1][j][k]  - phi[i-1][j][k] * (*u)[mc][i][j][k];
    
    mc = Comp::v();
    fluxy = phi[i][j][k] * (*u)[mc][i][j+1][k]  - phi[i][j-1][k] * (*u)[mc][i][j][k];

    mc = Comp::w();
    fluxz = phi[i][j][k] * (*u)[mc][i][j][k+1]  - phi[i][j][k-1] * (*u)[mc][i][j][k];

/* update clr */

    clr[i][j][k] = -fluxx * dtdx + phi[i][j][k];
   // std::cout<< clrn[i][j][k] <<std::ends;
}

  clr.bnd_update();
  clr.exchange_all();

  for_aijk(i,j,k)
    phi[i][j][k]=clr[i][j][k];

#endif

    phi.bnd_update();
    phi.exchange_all();

      std::cout<< dt <<" "<<dx <<"\n";
//    std::cout<< phi[80][1][1] <<"\n";

 /********************************************************+
 |              Final treatment                           |
 *+********************************************************/
//  for_aijk(i,j,k)
//    fext[i][j][k]=clr[i][j][k];

  return;

}
   



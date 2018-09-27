#include "vof.h"
#include <cmath>

/******************************************************************************/
void VOF::plic_advance() {

insert_bc(clr);
clr.exchange_all();

const real dt= time->dt();
real dx, dy, dz;
real dtdx, dtdy, dtdz;
real dx2, dx3, dy2, dy3, dz2, dz3;
real cx, cy, cz;
real fluxx, fluxy,fluxz;
Comp mc;

for_ijk(i,j,k){

    dx=clr.dxc(i);
    dy=clr.dyc(j);
    dz=clr.dzc(k);

    dtdx=dt/dx;
    dtdy=dt/dy;
    dtdz=dt/dz;

    dx2 =dx*dx;
    dx3 =dx2*dx;
    dy2 =dy*dy;
    dy3 =dy2*dy;
    dz2 =dz*dz;
    dz3 =dz2*dz;


/* Flux */
    mc = Comp::u();
    fluxx = clr[i][j][k] * (*u)[mc][i+1][j][k] * dtdx - clr[i-1][j][k] * (*u)[mc][i][j][k];
    
    mc = Comp::v();
    fluxy = clr[i][j][k] * (*u)[mc][i][j+1][k] * dtdy - clr[i][j-1][k] * (*u)[mc][i][j][k];

    mc = Comp::w();
    fluxz = clr[i][j][k] * (*u)[mc][i][j][k+1] * dtdz - clr[i][j][k-1] * (*u)[mc][i][j][k];

/* update clr */

    clrn[i][j][k] = -fluxx * dtdx + clr[i][j][k];

}

  insert_bc(clrn);
  clrn.exchange_all();

  for_aijk(i,j,k)
    clr[i][j][k]=clrn[i][j][k];

 /********************************************************+
 |              Final treatment                           |
 *+********************************************************/
  for_aijk(i,j,k)
    fext[i][j][k]=clr[i][j][k];

  return;

}
   



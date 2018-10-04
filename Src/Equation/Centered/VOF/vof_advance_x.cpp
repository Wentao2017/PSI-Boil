#include "vof.h"

/******************************************************************************/
void VOF::advance_x() {

  // calculate normal vector
  gradphic(phi);

  // advance in the x-direction

  Comp m = Comp::u();
  for_vmijk((*u),m,i,j,k){

    // upwind i-index
    int iup = i-1;
    if((*u)[m][i][j][k]<0.0) iup = i;

    // color function upwind
    real c = phi[iup][j][k];

    // calculate vn1, vn2, vn3: normal vector at face center
    real vn1, vn2, vn3;
    vn1 = nx[iup][j][k];
    vn2 = ny[iup][j][k];
    vn3 = nz[iup][j][k];

    // calculate g: CFL upwind
    real g;
    real dt=time->dt();
    g = ((*u)[m][i][j][k])*dt/phi.dxc(i-1);
    if((*u)[m][i][j][k]<0.0) g = ((*u)[m][i][j][k])*dt/phi.dxc(i);

   real absg = fabs(g);
   real vm1 = fabs(vn1);
   real vm2 = fabs(vn2);
   real vm3 = fabs(vn3)+boil::pico;
   real qa = 1.0/(vm1+vm2+vm3);
   vm1 *= qa;
   vm2 *= qa;
   vm3 *= qa;
   real alpha = calc_alpha(c, vm1, vm2, vm3);
   real ra = vm1 * (1.0 - absg);
   qa = 1.0/(1.0-ra);
   if (g*vn1 > 0) alpha = alpha -ra;
   vm1 = vm1 * absg;

   // calculate f: flux
    real f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g;

    // update color function store as stmp
    stmp[i-1][j][k] = phi[i-1][j][k]-f;
    stmp[i  ][j][k] = phi[i  ][j][k]+f;
  }

  // update phi
  phi = stmp;
  phi.bnd_update();
  phi.exchange();
}


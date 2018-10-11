#include "vof.h"

/******************************************************************************/
void VOF::advance_y() {

  // advance in the y-direction

  Comp m = Comp::v();
  for_vmijk((*u),m,i,j,k){

    // flux
    real f;

    // upwind j-index
    int jup = j-1;
    if((*u)[m][i][j][k]<0.0) jup = j;             

    // calculate g: CFL upwind
    real g;
    real dt=time->dt();
    g = ((*u)[m][i][j][k])*dt/phi.dyc(j-1);
    if((*u)[m][i][j][k]<0.0) g = ((*u)[m][i][j][k])*dt/phi.dyc(j);

    if (phi[i][jup][k]<boil::pico) {

      f = phi[i][jup][k] * g;
      //f = 0.0 * g;

    } else if(phi[i][jup][k]<1.0-boil::pico) {

      f = phi[i][jup][k] * g;
      //f = 1.0 * g;

    } else {

      // color function upwind
      real c = phi[i][jup][k];

      // calculate vn1, vn2, vn3: normal vector at face center
      real vn1, vn2, vn3;
      vn1 = -nx[i][jup][k];
      vn2 = -ny[i][jup][k];
      vn3 = -nz[i][jup][k];

      real absg = fabs(g);
      real vm1 = fabs(vn1);
      real vm2 = fabs(vn2);
      real vm3 = fabs(vn3)+boil::pico;
      real qa = 1.0/(vm1+vm2+vm3);
      vm1 *= qa;
      vm2 *= qa;
      vm3 *= qa;
      real alpha = calc_alpha(c, vm1, vm2, vm3);
      
      real ra = vm2 * (1.0 - absg);
      qa = 1.0/(1.0-ra);
      if (g*vn2 > 0) alpha = alpha -ra;
      vm2 = vm2 * absg;

      // calculate f: flux
      f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g;

    }

    // update color function store as stmp
    flx[m][i][j][k] = f * dV(i,jup,k);

  }

  //boil::plot->plot(flux_x, "flux_x", time->current_step());
 
  return;
}


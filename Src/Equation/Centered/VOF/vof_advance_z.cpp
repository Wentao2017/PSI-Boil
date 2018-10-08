#include "vof.h"

/******************************************************************************/
void VOF::advance_z() {

  // calculate normal vector
  gradphic(phi);
  
//boil::plot->plot(nx, ny, nz, "nx-ny-nz", time->current_step());
#if 0
  std::cout<<"nx50 "<<nx[50][1][1]<<" "<<"nx51 "<<nx[51][1][1]<<"\n";
  std::cout<<"ny50 "<<ny[50][1][1]<<" "<<"ny51 "<<ny[51][1][1]<<"\n";
  std::cout<<"nz50 "<<nz[50][1][1]<<" "<<"nz51 "<<nz[51][1][1]<<"\n";
  std::cout<<"nx75 "<<nx[75][1][1]<<" "<<"nx76 "<<nx[76][1][1]<<"\n";
  std::cout<<"ny75 "<<ny[75][1][1]<<" "<<"ny76 "<<ny[76][1][1]<<"\n";
  std::cout<<"nz75 "<<nz[75][1][1]<<" "<<"nz76 "<<nz[76][1][1]<<"\n";
#endif
 
  // advance in the z-direction
  for_ijk(i,j,k){
    stmp[i][j][k]=phi[i][j][k] * dV(i,j,k);
  }

  Comp m = Comp::w();
  for_vmijk((*u),m,i,j,k){

    // flux
    real f;

    // upwind j-index
    int kup = k-1;
    if((*u)[m][i][j][k]<0.0) kup = k;             

    // calculate g: CFL upwind
    real g;
    real dt=time->dt();
    g = ((*u)[m][i][j][k])*dt/phi.dzc(k-1);
    if((*u)[m][i][j][k]<0.0) g = ((*u)[m][i][j][k])*dt/phi.dzc(k);

    if (approx(phi[i][j][kup], 0.0)) {

      f = 0.0 * g;

    } else if(approx(phi[i][j][kup],1.0)) {

      f = 1.0 * g;

    } else {

      // color function upwind
      real c = phi[i][j][kup];

      // calculate vn1, vn2, vn3: normal vector at face center
      real vn1, vn2, vn3;
      vn1 = -nx[i][j][kup];
      vn2 = -ny[i][j][kup];
      vn3 = -nz[i][j][kup];

      real absg = fabs(g);
      real vm1 = fabs(vn1);
      real vm2 = fabs(vn2);
      real vm3 = fabs(vn3)+boil::pico;
      real qa = 1.0/(vm1+vm2+vm3);
      vm1 *= qa;
      vm2 *= qa;
      vm3 *= qa;
      real alpha = calc_alpha(c, vm1, vm2, vm3);
      
      real ra = vm3 * (1.0 - absg);
      qa = 1.0/(1.0-ra);
      if (g*vn3 > 0) alpha = alpha -ra;
      vm3 = vm3 * absg;

      // calculate f: flux
      f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g;

    }

    // update color function store as stmp
    stmp[i][j][k-1] = stmp[i][j][k-1]-f * dV(i,j,kup);
    stmp[i  ][j][k] = stmp[i  ][j][k]+f * dV(i,j,kup);

  }

  // update phi
  for_ijk(i,j,k){
    real phi_tmp = stmp[i][j][k] / dV(i,j,k);
    phi[i][j][k] = std::min(1.0,std::max(0.0,phi_tmp));
  }
  phi.bnd_update();
  phi.exchange_all();

  //boil::plot->plot(flux_x, "flux_x", time->current_step());
 
  return;
}


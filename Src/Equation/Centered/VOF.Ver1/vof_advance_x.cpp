#include "vof.h"

/******************************************************************************/
void VOF::advance_x() {

  // calculate normal vector
  gradphic(stmp);
  
//boil::plot->plot(nx, ny, nz, "nx-ny-nz", time->current_step());
#if 0
  std::cout<<"nx50 "<<nx[50][1][1]<<" "<<"nx51 "<<nx[51][1][1]<<"\n";
  std::cout<<"ny50 "<<ny[50][1][1]<<" "<<"ny51 "<<ny[51][1][1]<<"\n";
  std::cout<<"nz50 "<<nz[50][1][1]<<" "<<"nz51 "<<nz[51][1][1]<<"\n";
  std::cout<<"nx75 "<<nx[75][1][1]<<" "<<"nx76 "<<nx[76][1][1]<<"\n";
  std::cout<<"ny75 "<<ny[75][1][1]<<" "<<"ny76 "<<ny[76][1][1]<<"\n";
  std::cout<<"nz75 "<<nz[75][1][1]<<" "<<"nz76 "<<nz[76][1][1]<<"\n";
#endif
 
  // advance in the x-direction
  for_ijk(i,j,k){
    stmp[i][j][k]=phi[i][j][k] * dV(i,j,k);
  }

  Comp m = Comp::u();
//  std::cout<<"Debug 0 "<<"\n";
  for_vmijk((*u),m,i,j,k){

    // flux
    real f;

    // upwind i-index
    int iup = i-1;
    if((*u)[m][i][j][k]<0.0) iup = i;             

    // calculate g: CFL upwind
    real g;
    real dt=time->dt();
    g = ((*u)[m][i][j][k])*dt/phi.dxc(i-1);
    if((*u)[m][i][j][k]<0.0) g = ((*u)[m][i][j][k])*dt/phi.dxc(i);

    //if (approx(phi[iup][j][k], 0.0)) {
    if (phi[iup][j][k] < boil::pico) {

      f = phi[iup][j][k] * g;
      //f = 0.0 * g;

    //} else if(approx(phi[iup][j][k],1.0)) {
    } else if(phi[iup][j][k]>1.0-boil::pico) {

      f = phi[iup][j][k] * g;
      //f = 1.0 * g;

    } else {

//      std::cout<<"iup= "<<iup<<" j= "<<j<<" k= "<<k<<"\n";
      // color function upwind
      real c = phi[iup][j][k];

      // calculate vn1, vn2, vn3: normal vector at face center
      real vn1, vn2, vn3;
      vn1 = -nx[iup][j][k];
      vn2 = -ny[iup][j][k];
      vn3 = -nz[iup][j][k];

  //    std::cout<<"n= "<<vn1<<" "<<vn2<<" "<<vn3<<"\n";

      real absg = fabs(g);
      real vm1 = fabs(vn1);
      real vm2 = fabs(vn2);
      real vm3 = fabs(vn3)+boil::pico;
      real qa = 1.0/(vm1+vm2+vm3);
      vm1 *= qa;
      vm2 *= qa;
      vm3 *= qa;
    //  std::cout<<"Test 0 "<<"\n";
    //  std::cout<<"i= "<<i-1<<" j= "<<j-1<<" k= "<<k-1<<"\n";
      real alpha = calc_alpha(c, vm1, vm2, vm3);
      
      real ra = vm1 * (1.0 - absg);
      qa = 1.0/(1.0-ra);
      if (g*vn1 > 0) alpha = alpha -ra;
      vm1 = vm1 * absg;

      // calculate f: flux
      f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g;

    }

    // update color function store as stmp
    stmp[i-1][j][k] = stmp[i-1][j][k]-f * dV(iup,j,k);
    stmp[i  ][j][k] = stmp[i  ][j][k]+f * dV(iup,j,k);

  }

  // update phi
  for_ijk(i,j,k){
    phi[i][j][k] = stmp[i][j][k] / dV(i,j,k);
  }
  phi.bnd_update();
  phi.exchange_all();

  // update stmp
  for_ijk(i,j,k){
    stmp[i][j][k] = std::min(1.0,std::max(0.0,phi[i][j][k]));
  }
  stmp.bnd_update();
  stmp.exchange_all();

  //boil::plot->plot(flux_x, "flux_x", time->current_step());
 
  return;
}


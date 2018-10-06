#include "vof.h"

/******************************************************************************/
void VOF::advance_x() {

  // calculate normal vector
  gradphic(phi);
  
  boil::plot->plot(nx, ny, nz, "nx-ny-nz", time->current_step());
  std::cout<<"nx50 "<<nx[50][1][1]<<" "<<"nx51 "<<nx[51][1][1]<<"\n";
  std::cout<<"ny50 "<<ny[50][1][1]<<" "<<"ny51 "<<ny[51][1][1]<<"\n";
  std::cout<<"nz50 "<<nz[50][1][1]<<" "<<"nz51 "<<nz[51][1][1]<<"\n";
  std::cout<<"nx75 "<<nx[75][1][1]<<" "<<"nx76 "<<nx[76][1][1]<<"\n";
  std::cout<<"ny75 "<<ny[75][1][1]<<" "<<"ny76 "<<ny[76][1][1]<<"\n";
  std::cout<<"nz75 "<<nz[75][1][1]<<" "<<"nz76 "<<nz[76][1][1]<<"\n";
 
  // advance in the x-direction

  stmp=phi;

  Comp m = Comp::u();
  for_vmijk((*u),m,i,j,k){

    //std::cout<<i<<" "<<j<<" "<<k<<"\n";

    // flux
    real f;

    // calculate g: CFL upwind
    real g;
    real dt=time->dt();
    g = ((*u)[m][i][j][k])*dt/phi.dxc(i-1);
    if((*u)[m][i][j][k]<0.0) g = ((*u)[m][i][j][k])*dt/phi.dxc(i);

    if (approx(phi[i-1][j][k], 0.0) && approx(phi[i][j][k], 0.0)) {
      f = 0.0;
    } else if(approx(phi[i-1][j][k],1.0)&&approx(phi[i][j][k],1.0)) {
      f = g;
    } else {

      // upwind i-index
      int iup = i-1;
      if((*u)[m][i][j][k]<0.0) iup = i;

      // color function upwind
      real c = phi[iup][j][k];

      // calculate vn1, vn2, vn3: normal vector at face center
      real vn1, vn2, vn3;
      vn1 = -nx[iup][j][k];
      vn2 = -ny[iup][j][k];
      vn3 = -nz[iup][j][k];

      real absg = fabs(g);
      real vm1 = fabs(vn1);
      real vm2 = fabs(vn2);
      real vm3 = fabs(vn3)+boil::pico;
      real qa = 1.0/(vm1+vm2+vm3);
      vm1 *= qa;
      vm2 *= qa;
      vm3 *= qa;
      
      std::cout<<"i= "<<i-1<<" j= "<<j<<" k= "<<k<<"\n"; 
      std::cout<<"phi[50] "<<phi[50][1][1]<<" phi[51] "<<phi[51][1][1] <<"\n";
 
      real alpha = calc_alpha(c, vm1, vm2, vm3);
      
      if(i==52&&j==1&&k==1){
        std::cout<<"alpha[52] "<<alpha<<"\n";
      }

      real ra = vm1 * (1.0 - absg);
      qa = 1.0/(1.0-ra);
      if (g*vn1 > 0) alpha = alpha -ra;
      vm1 = vm1 * absg;

      if(i==52&&j==1&&k==1){
        std::cout<<"alpha[52] "<<alpha * qa<<"\n";
      }

       
       // calculate f: flux
      f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g;
    }

    // update color function store as stmp
    stmp[i-1][j][k] = stmp[i-1][j][k]-f;
    stmp[i  ][j][k] = stmp[i  ][j][k]+f;

#if 0
    if((i==75||i==76)&&j==1&&k==1){
      std::cout<<"\n";
      //std::cout<<"iup= "<<iup<<" c= "<<c<<" vn1= "<<vn1<<" vn2= "<<vn2<<" vn3= "<<vn3<<"\n";
      //std::cout<<"f= "<<f<<" g= "<<g<<" vm1= "<<vm1<<" vm2= "<<vm2<<" vm3= "<<vm3<<"\n";
      std::cout<<"f= "<<f<<"\n";
      //std::cout<<"alpha= "<<alpha<<" calc_v= "<<f/g<<"\n";
      std::cout<<"stmp[i-1]= "<<stmp[i-1][j][k]<<" "<<phi[i-1][j][k]<<"\n";
      std::cout<<"stmp[i  ]= "<<stmp[i  ][j][k]<<" "<<phi[i  ][j][k]<<"\n";
      std::cout<<"\n";
    }
#endif

  }

#if 0
  std::cout<<"stmp[75]= "<<stmp[75][1][1]<<" "<<stmp[76][1][1]<<"\n";
#endif

  // update phi
  phi = stmp;
  phi.bnd_update();
  phi.exchange_all();

  //boil::plot->plot(flux_x, "flux_x", time->current_step());
 
  return;
}


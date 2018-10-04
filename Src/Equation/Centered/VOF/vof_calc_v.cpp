#include "vof.h"
#include <algorithm>

#if 0
void VOF::calc_v(){

  for_ijk(i,j,k){
 
    a[i][j][k]   = boil::minr(alpha[i][j][k], 1.0-alpha[i][j][k]);
    vv[i][j][k]  = 0.0;

    if(a[i][j][k] > 0){
      
      vm1[i][j][k]  = boil::minr(vma[i][j][k], vmb[i][j][k], vmc[i][j][k]);
      vm3[i][j][k]  = boil::maxr(vma[i][j][k], vmb[i][j][k], vmc[i][j][k]);
      vm2[i][j][k]  = fabs(1.0 - vm3[i][j][k] - vm1[i][j][k]);
      vm12[i][j][k] = vm1[i][j][k] + vm2[i][j][k];
    
      if(a[i][j][k] < vm1[i][j][k]){
        
        vv[i][j][k] = pow(a[i][j][k], 3.0) / (6.0 * vm1[i][j][k] * vm2[i][j][k] * vm3[i][j][k]);
   
      }
      else if(a[i][j][k] < vm2[i][j][k]){
 
        vv[i][j][k] = a[i][j][k] * (a[i][j][k] - vm1[i][j][k]) / (2.0 * vm2[i][j][k] * vm3[i][j][k]) + \
                      pow(vm1[i][j][k], 2.0) / (6.0 * vm2[i][j][k] * vm3[i][j][k] + boil::pico);

      }
      else if(a[i][j][k] < boil::minr(vm12[i][j][k], vm3[i][j][k])){

        vv[i][j][k] = (pow(a[i][j][k], 2.0) * (3.0 * vm12[i][j][k] - a[i][j][k]) + pow(vm1[i][j][k], 2.0) * \
                      (vm1[i][j][k] - 3.0 * a[i][j][k]) + pow(vm2[i][j][k], 2.0) * (vm2[i][j][k] -3.0 * a[i][j][k])) / \
                      (6.0 * vm1[i][j][k] * vm2[i][j][k] * vm3[i][j][k]);

      }
      else if(vm3[i][j][k] < vm12[i][j][k]){

        vv[i][j][k] = (pow(a[i][j][k], 2.0) * (3.0 - 2.0 * a[i][j][k]) + pow(vm1[i][j][k], 2.0) * (vm1[i][j][k] - 3.0 * a[i][j][k]) + \
                      pow(vm2[i][j][k], 2.0) * (vm2[i][j][k] - 3.0 * a[i][j][k]) + pow(vm3[i][j][k], 2.0) * (vm3[i][j][k] - 3.0 * \
                      a[i][j][k])) / (6.0 * vm1[i][j][k] * vm2[i][j][k] * vm3[i][j][k]);
    
      }
      else{

        vv[i][j][k] = (a[i][j][k] - 0.5 * vm12[i][j][k]) / vm3[i][j][k];

      }
    }
    
    if(alpha[i][j][k] > 0.5){
    
      vv[i][j][k] = 1.0 - vv[i][j][k];
    
    }
  }
  vv.bnd_update();
  vv.exchange_all();

  return;

}
#endif

real VOF::calc_v(real alpha, real vma, real vmb, real vmc){
  real v;
  return v;
}

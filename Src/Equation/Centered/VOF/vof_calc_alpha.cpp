#include "vof.h"
//#include <algorithm>

#if 0
void VOF::calc_alpha(){

    for_ijk(i,j,k){

      w[i][j][k]    = boil::minr(phi[i][j][k], 1.0-phi[i][j][k]);
      vm1[i][j][k]  = boil::minr(vma[i][j][k], vmb[i][j][k], vmc[i][j][k]);     
      vm3[i][j][k]  = boil::maxr(vma[i][j][k], vmb[i][j][k], vmc[i][j][k]);
      vm2[i][j][k]  = fabs(1.0 - vm3[i][j][k] - vm1[i][j][k]);
      vm12[i][j][k] = vm1[i][j][k] + vm2[i][j][k];
    
      v1[i][j][k] = pow(vm1[i][j][k],2.0) / (6.0 * vm2[i][j][k] * vm3[i][j][k] + boil::pico);
    
    }

    std::ofstream fout0;
    fout0.open("w.txt");
      for_aijk(i,j,k)
      fout0<<"w "<<i<<" "<<j<<" "<<k<<" "<<w[i][j][k] <<"\n";
    fout0.close();
    
    std::ofstream fout1;
    fout1.open("vm1.txt");
      for_aijk(i,j,k)
      fout1<<"vm1 "<<i<<" "<<j<<" "<<k<<" "<<vm1[i][j][k] <<"\n";
    fout1.close();

    std::ofstream fout2;
    fout2.open("vm2.txt");
      for_aijk(i,j,k)
      fout2<<"vm2 "<<i<<" "<<j<<" "<<k<<" "<<vm2[i][j][k] <<"\n";
    fout2.close();

    std::ofstream fout3;
    fout3.open("vm3.txt");
      for_aijk(i,j,k)
      fout3<<"vm3 "<<i<<" "<<j<<" "<<k<<" "<<vm3[i][j][k] <<"\n";
    fout3.close();
 
    std::ofstream fout4;
    fout4.open("vm12.txt");
      for_aijk(i,j,k)
      fout4<<"vm12 "<<i<<" "<<j<<" "<<k<<" "<<vm12[i][j][k] <<"\n";
    fout4.close();
  
    std::ofstream fout5;
    fout5.open("v1.txt");
      for_aijk(i,j,k)
      fout5<<"v1 "<<i<<" "<<j<<" "<<k<<" "<<v1[i][j][k] <<"\n";
    fout5.close();

    std::ofstream fout6;
    fout6.open("phi.txt");
      for_aijk(i,j,k)
      fout6<<"phi "<<i<<" "<<j<<" "<<k<<" "<<phi[i][j][k] <<"\n";
    fout6.close();



    boil::plot->plot(w, "w");
    boil::plot->plot(vm1, "vm1");
    boil::plot->plot(vm2, "vm2");
    boil::plot->plot(vm3, "vm3");
    boil::plot->plot(vm12, "vm12");
    boil::plot->plot(v1, "v1");
    boil::plot->plot(phi, "phi");

    return;
}
#endif



#if 1
void VOF::calc_alpha(){


  for_ijk(i,j,k){

    w[i][j][k]    = boil::minr(phi[i][j][k], 1.0-phi[i][j][k]);
    vm1[i][j][k]  = boil::minr(vma[i][j][k], vmb[i][j][k], vmc[i][j][k]);     
    vm3[i][j][k]  = boil::maxr(vma[i][j][k], vmb[i][j][k], vmc[i][j][k]);
    vm2[i][j][k]  = fabs(1.0 - vm3[i][j][k] - vm1[i][j][k]);
    vm12[i][j][k] = vm1[i][j][k] + vm2[i][j][k];
    
    v1[i][j][k] = pow(vm1[i][j][k],2.0) / (6.0 * vm2[i][j][k] * vm3[i][j][k] + boil::pico);
    
    if(w[i][j][k] < v1[i][j][k]){

      alpha[i][j][k] = pow(6.0 * vm1[i][j][k] * vm2[i][j][k] * vm3[i][j][k] * w[i][j][k], 1.0/3);

    }
    else if (w[i][j][k] < v1[i][j][k] + (vm2[i][j][k] - vm1[i][j][k]) / (2.0 * vm3[i][j][k])){
    
      alpha[i][j][k] = 0.5 * (vm1[i][j][k] + sqrt(pow(vm1[i][j][k],2) + 8.0 * vm2[i][j][k] * vm3[i][j][k] * (w[i][j][k] - v1[i][j][k])));
 
    }
    else{
    
      alpha[i][j][k] = 0.0;
      if (vm3[i][j][k] < vm12[i][j][k]){

        v3[i][j][k] = (pow(vm3[i][j][k],2.0) * (3.0 * vm12[i][j][k] - vm3[i][j][k]) + pow(vm1[i][j][k],2.0) * (vm1[i][j][k] - 3.0 * vm3[i][j][k]) + \
                      pow(vm2[i][j][k],2.0) * (vm2[i][j][k] - 3.0 * vm3[i][j][k])) / (6.0 * vm1[i][j][k] * vm2[i][j][k] * vm3[i][j][k]);

      }
      else{
   
        v3[i][j][k] = 0.5 * vm12[i][j][k] / vm3[i][j][k];
        if(v3[i][j][k] <= w[i][j][k]){

          alpha[i][j][k] = vm3[i][j][k] * w[i][j][k] + 0.5 * vm12[i][j][k];

        }            
      }
      if (alpha[i][j][k] == 0.0){
        
        if(w[i][j][k] < v3[i][j][k]){

          a2[i][j][k] = -3.0 * vm12[i][j][k];
          a1[i][j][k] = 3.0 * (pow(vm1[i][j][k], 2.0) + pow(vm2[i][j][k], 2.0));               
          a0[i][j][k] = -(pow(vm1[i][j][k], 3.0) + pow(vm2[i][j][k], 3.0) - 6.0 * vm1[i][j][k] * vm2[i][j][k] * vm3[i][j][k] * w[i][j][k]);
        
        }
        else{

          a2[i][j][k] = -1.5;
          a1[i][j][k] = 1.5 * (pow(vm1[i][j][k], 2.0) + pow(vm2[i][j][k], 2.0) + pow(vm3[i][j][k], 2.0));
          a0[i][j][k] = -0.5 * (pow(vm1[i][j][k], 3.0) + pow(vm2[i][j][k], 3.0) + pow(vm3[i][j][k], 3.0) - 6.0 * vm1[i][j][k] * vm2[i][j][k] * \
                       vm3[i][j][k] * w[i][j][k]);
        }
        q0[i][j][k]  = (1.0/6) * (a1[i][j][k] * a2[i][j][k] - 3.0 * a0[i][j][k]) - pow(a2[i][j][k], 3.0) * (1.0/27);
        sp[i][j][k]  = sqrt(-1.0/3 * a1[i][j][k] + 1.0/9 * pow(a2[i][j][k], 2.0));
        th[i][j][k]  = 1.0/3 * acos(q0[i][j][k] / (pow(sp[i][j][k], 3.0)));
        alpha[i][j][k] = 2.0 * sp[i][j][k] * cos(th[i][j][k] + (4.0/3 * (boil::pi))) - (1.0/3) * a2[i][j][k];
    
      }
    } 
    if(phi[i][j][k] > 0.5){
  
      alpha[i][j][k] = 1.0 - alpha[i][j][k];

    }
  }
  
  alpha.bnd_update();
  alpha.exchange_all();
  
  #if 1
  std::ofstream fout0;
  fout0.open("alpha.txt");
  for_aijk(i,j,k)
    fout0<<"alpha "<<i<<" "<<j<<" "<<k<<" "<<alpha[i][j][k] <<"\n";
  fout0.close();
  #endif

  return;

}

real VOF::calc_alpha(real & v, real & vma, real & vmb, real & vmc){
  real w = boil::minr(v, 1.0-v);
  real alpha;
  return alpha;
}
#endif

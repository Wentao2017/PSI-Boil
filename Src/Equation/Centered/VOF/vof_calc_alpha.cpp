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



#if 0
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
#endif

real VOF::calc_alpha(real & v, real & vma, real & vmb, real & vmc){

  real w, v1, v3, vm1, vm2, vm3, vm12; 
  real a0, a1, a2, q0, sp, th;
  real alpha;

  std::cout<<"calc_alpha:point00 \n";

  w     = boil::minr(v, 1.0-v);  
  vm1   = boil::minr(vma, vmb, vmc);
  vm3   = boil::maxr(vma, vmb, vmc);
  vm2   = fabs(1.0 - vm3 - vm1);
  vm12  = vm1 + vm2;

  v1    = pow(vm1, 2.0) / (6.0 * vm2 * vm3 + boil::pico);

  std::cout<<"calc_alpha:point10 \n";

  if (w < v1){
    alpha = pow(6.0 * vm1 * vm2 * vm3 * w, 1.0/3.0);
  } else if(w < v1 + (vm2-vm1) / (2.0 * vm3)){
    alpha = 0.5 * (vm1 + sqrt(pow(vm1, 2) + 8.0 * vm2 * vm3 * (w - v1)));
  } else {
    std::cout<<"calc_alpha:point15 \n";
    alpha = 0.0;
    if (vm3 < vm12){
      std::cout<<"calc_alpha:point16 \n";
      v3 = (pow(vm3, 2.0) * (3.0 * vm12 - vm3) + pow(vm1, 2.0) * (vm1 - 3.0 *
           vm3) + pow(vm2, 2.0) * (vm2 - 3.0 * vm3)) / (6.0 * vm1 * vm2 * vm3);
    } else{
      std::cout<<"calc_alpha:point17 \n";
      v3 = 0.5 * vm12 / vm3;
      if(v3 <= w){
        alpha = vm3 * w + 0.5 * vm12;
      }
    }
    std::cout<<"calc_alpha:point18 \n";
    if (alpha == 0.0){
      std::cout<<"calc_alpha:point19 \n";
      if(w < v3){
        std::cout<<"calc_alpha:point20 \n";
        a2 = -3.0 * vm12;
        a1 = 3.0 * (pow(vm1, 2.0) + pow(vm2, 2.0));               
        a0 = -(pow(vm1, 3.0) + pow(vm2, 3.0) - 6.0 * vm1 * vm2 * vm3 * w);
      } else {
        std::cout<<"calc_alpha:point21 \n";
        a2 = -1.5;
        a1 = 1.5 * (pow(vm1, 2.0) + pow(vm2, 2.0) + pow(vm3, 2.0));
        a0 = -0.5 * (pow(vm1, 3.0) + pow(vm2, 3.0) + pow(vm3, 3.0) - 6.0 *
             vm1 * vm2 * vm3 * w);
      }
      std::cout<<"calc_alpha:point22 \n";
      q0     = (1.0/6) * (a1 * a2 - 3.0 * a0) - pow(a2, 3.0) * (1.0/27.0);
      std::cout<<"calc_alpha:point23 \n";
      std::cout<<"inside sp "<<-1.0/3.0 * a1 + 1.0/9.0 * pow(a2, 2.0)<<"\n";
      sp    = sqrt(-1.0/3.0 * a1 + 1.0/9.0 * pow(a2, 2.0));
      std::cout<<"calc_alpha:point24 \n";
      th    = 1.0/3.0 * acos(q0 / (pow(sp, 3.0)));
      std::cout<<"calc_alpha:point25 \n";
      alpha = 2.0 * sp * cos(th + (4.0/3.0 * (boil::pi))) - (1.0/3.0) * a2;
      std::cout<<"calc_alpha:point29 \n";
    }
  }
 
  //std::cout<<"calc_alpha:point30 \n";

  if(v > 0.5){
      alpha = 1.0 - alpha;
  }
  
  return alpha;
}


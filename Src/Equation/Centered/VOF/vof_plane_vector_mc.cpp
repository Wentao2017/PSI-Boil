#include "vof.h"
#include <fstream>
void VOF::plane_vector_mc() {
/***************************************************************************//**
*  \ vma +vmb+vmc =1, vma [0,1], vmb [0,1], vmc[0,1]
*  *******************************************************************************/ 

  real a, b, c;
  for_ijk(i,j,k){
    
    a = fabs(nx[i][j][k]) / (fabs(nx[i][j][k]) + fabs(ny[i][j][k]) + fabs (nz[i][j][k]) + boil::pico);
    b = fabs(ny[i][j][k]) / (fabs(nx[i][j][k]) + fabs(ny[i][j][k]) + fabs (nz[i][j][k]) + boil::pico);
    c = fabs(nz[i][j][k]) / (fabs(nx[i][j][k]) + fabs(ny[i][j][k]) + fabs (nz[i][j][k]) + boil::pico);
    
    vma[i][j][k] = a;
    vmb[i][j][k] = b;
    vmc[i][j][k] = c;
  }
    

  vma.exchange_all();
  vmb.exchange_all();
  vmc.exchange_all();

#if 0
  std::ofstream fout0;
  fout0.open("vma.txt");
  for_aijk(i,j,k)
    fout0<<"vma "<<i<<" "<<j<<" "<<k<<" "<<vma[i][j][k] <<"\n";
  fout0.close();

  std::ofstream fout1;
  fout1.open("vmb.txt");
  for_aijk(i,j,k)
    fout1<<"vmb "<<i<<" "<<j<<" "<<k<<" "<<vmb[i][j][k] <<"\n";
  fout1.close();

  std::ofstream fout2;
  fout2.open("vmc.txt");
  for_aijk(i,j,k)
    fout2<<"vmc "<<i<<" "<<j<<" "<<k<<" "<<vmc[i][j][k] <<"\n";
  fout2.close();
#endif


  return;

}

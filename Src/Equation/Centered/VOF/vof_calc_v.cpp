#include "vof.h"
#include <algorithm>

void VOF::calc_v(){

  for_ijk(i,j,k){
 
    a[i][j][k]   = boil::minr(alpha[i][j][k], 1.0-alpha[i][j][k]);
    vv[i][j][k]  = 0.0;

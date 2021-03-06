#include "vof.h"
#include <iomanip>

/******************************************************************************/
void VOF::set_iflag() {
/***************************************************************************//**
*  /brief flagging
*           output   :iflag
*
*   iflag: -3 | -2 | -1 |  0 |  0 |  1 |  2 |  3 | ..   
*          ---+----+----+----+----+----+----+----+---
*   vapor                   |                        liquid
*                      free-surface
*
*   solid region: iflag = -1001
*******************************************************************************/
  const int ifmax =  nlayer+4;
  const int ifmin = -nlayer-4;

  /*-----------+
  |  set flag  |
  +-----------*/
  for_aijk(i,j,k) {
    if(phi[i][j][k]<phisurf){
      iflag[i][j][k]=ifmin;
    } else {
      iflag[i][j][k]=ifmax;
    }
  }

  for_ijk(i,j,k) {
    if(dom->ibody().off(i,j,k))
      iflag[i][j][k]=-1001;
  }
  insert_bc_flag(iflag, false);
  iflag.exchange();

  /* next to free-surface (NFCell) */
  /* i-direction */
  for(int i=0; i<ni()-1; i++){
    for_jk(j,k){
      if ((phi[i][j][k]-phisurf)*(phi[i+1][j][k]-phisurf)<=0.0) {
        if(iflag[i  ][j][k]>-1000) iflag[i  ][j][k]=0;
        if(iflag[i+1][j][k]>-1000) iflag[i+1][j][k]=0;
      }
    }
  }
  /* j-direction */
  for(int j=0; j<nj()-1; j++){
    for_ik(i,k){
      if((phi[i][j][k]-phisurf)*(phi[i][j+1][k]-phisurf)<=0.0){
        if(iflag[i][j  ][k]>-1000) iflag[i][j  ][k]=0;
        if(iflag[i][j+1][k]>-1000) iflag[i][j+1][k]=0;
      }
    }
  }
  /* k-direction */
  for(int k=0; k<nk()-1; k++){
    for_ij(i,j){
      if((phi[i][j][k]-phisurf)*(phi[i][j][k+1]-phisurf)<=0.0){
        if(iflag[i][j][k  ]>-1000) iflag[i][j][k  ]=0;
        if(iflag[i][j][k+1]>-1000) iflag[i][j][k+1]=0;
      }
    }
  }
  insert_bc_flag(iflag, false);
  iflag.exchange();

  for(int layer=1; layer<=nlayer; layer++){
    /* i-direction */
    for(int i=0; i<ni()-1; i++){
      for_jk(j,k){
        if(int(abs(iflag[i  ][j][k]))==ifmax &&
           int(abs(iflag[i+1][j][k]))==(layer-1)){
           if(iflag[i  ][j][k]>-1000)
             iflag[i  ][j][k]=layer*int(copysign(1.0,iflag[i  ][j][k]));
        } else if(int(abs(iflag[i+1][j][k]))==ifmax &&
                  int(abs(iflag[i  ][j][k]))==(layer-1)){
           if(iflag[i+1][j][k]>-1000)
             iflag[i+1][j][k]=layer*int(copysign(1.0,iflag[i+1][j][k]));
        }
      }
    }
    /* j-direction */
    for(int j=0; j<nj()-1; j++){
      for_ik(i,k){
        if(int(abs(iflag[i][j  ][k]))==ifmax &&
           int(abs(iflag[i][j+1][k]))==(layer-1)){
           if(iflag[i][j  ][k]>-1000)
             iflag[i][j  ][k]=layer*int(copysign(1.0,iflag[i  ][j][k]));
        } else if(int(abs(iflag[i][j+1][k]))==ifmax &&
                  int(abs(iflag[i][j  ][k]))==(layer-1)){
           if(iflag[i][j+1][k]>-1000)
             iflag[i][j+1][k]=layer*int(copysign(1.0,iflag[i][j+1][k]));
        }
      }
    }
    /* k-direction */
    for(int k=0; k<nk()-1; k++){
      for_ij(i,j){
        if(int(abs(iflag[i][j][k  ]))==ifmax &&
           int(abs(iflag[i][j][k+1]))==(layer-1)){
           if(iflag[i][j][k  ]>-1000)
             iflag[i][j][k  ]=layer*int(copysign(1.0,iflag[i  ][j][k]));
        } else if(int(abs(iflag[i][j][k+1]))==ifmax &&
                  int(abs(iflag[i][j][k  ]))==(layer-1)){
           if(iflag[i][j][k+1]>-1000)
             iflag[i][j][k+1]=layer*int(copysign(1.0,iflag[i][j][k+1]));
        }
      }
    }
    insert_bc_flag(iflag, false);
    iflag.exchange();
  }

#if 0
  boil::plot->plot(phi, iflag, "phi-iflag", time->current_step());
  //exit(0);
#endif

  return;
}

/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_set_iflag.cpp,v 1.1 2015/05/06 08:00:09 sato Exp $'/
+-----------------------------------------------------------------------------*/

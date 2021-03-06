#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::norm_cc(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Centered Columns method in E.Aulisa,JCP,225(2007),2301-2319
*         Resluts: nx, ny, nz
*******************************************************************************/

  for_ijk(i,j,k) {
    real mxX, myX, mzX;
    mxX = copysign(1.0,+(sca[i+1][j][k]-sca[i-1][j][k]));
    myX = 0.5 * ( (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k])
                - (sca[i+1][j-1][k]+sca[i][j-1][k]+sca[i-1][j-1][k])); 
    mzX = 0.5 * ( (sca[i+1][j][k+1]+sca[i][j][k+1]+sca[i-1][j][k+1])
                - (sca[i+1][j][k-1]+sca[i][j][k-1]+sca[i-1][j][k-1])); 
    normalize(mxX,myX,mzX);

    real mxY, myY, mzY;
    mxY = 0.5 * ( (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
                - (sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k])); 
    myY = copysign(1.0,+(sca[i][j+1][k]-sca[i][j-1][k]));
    mzY = 0.5 * ( (sca[i][j-1][k+1]+sca[i][j][k+1]+sca[i][j+1][k+1])
                - (sca[i][j-1][k-1]+sca[i][j][k-1]+sca[i][j+1][k-1]));
    normalize(mxY,myY,mzY);

    real mxZ, myZ, mzZ;
    mxZ = 0.5 * ( (sca[i+1][j][k-1]+sca[i+1][j][k]+sca[i+1][j][k+1])
                - (sca[i-1][j][k-1]+sca[i-1][j][k]+sca[i-1][j][k+1])); 
    myZ = 0.5 * ( (sca[i][j+1][k-1]+sca[i][j+1][k]+sca[i][j+1][k+1])
                - (sca[i][j-1][k-1]+sca[i][j-1][k]+sca[i][j-1][k+1])); 
    mzZ = copysign(1.0,+(sca[i][j][k+1]-sca[i][j][k-1]));
    normalize(mxZ,myZ,mzZ);

    if (fabs(mxX)<fabs(myY)) {
      if (fabs(myY)<fabs(mzZ)) {
        nx[i][j][k]=mxZ;
        ny[i][j][k]=myZ;
        nz[i][j][k]=mzZ;
      } else {
        nx[i][j][k]=mxY;
        ny[i][j][k]=myY;
        nz[i][j][k]=mzY;
      }
    } else {
      if (fabs(mxX)<fabs(mzZ)) {
        nx[i][j][k]=mxZ;
        ny[i][j][k]=myZ;
        nz[i][j][k]=mzZ;
      } else {
        nx[i][j][k]=mxX;
        ny[i][j][k]=myX;
        nz[i][j][k]=mzX;
      }
    }
  }

  /* normal vector at adjacent cells next to wall, symmetric and IB */
  //insert_bc_gradphic(sca); 

  /* normal vector on boundary plane */
  insert_bc_norm_cc(sca);

  /* normalize */
  //for_avijk(sca,i,j,k) {
  //  normalize(nx[i][j][k],ny[i][j][k],nz[i][j][k]);
  //}

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //exit(0);


  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_gradphic.cpp,v 1.7 2015/05/05 15:11:07 sato Exp $'/
+-----------------------------------------------------------------------------*/

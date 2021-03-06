#include "vof.h"

/******************************************************************************/
void VOF::totalvol() {

   /*---------+
   | method 1 |
   +---------*/
   real phisum = 0.0;

   for_ijk(i,j,k)
     phisum += phi[i][j][k] * dV(i,j,k);

   boil::cart.sum_real(&phisum);


   std::cout.setf(std::ios_base::scientific);
   boil::oout << "totalvol:time,volume,phisum= " 
              << time->current_time()
              <<" "<< phisum << boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);

   return;
}
/*-----------------------------------------------------------------------------+
 '$Id: colorcip_totalvol.cpp,v 1.3 2009/11/12 12:15:48 sato Exp $'/
+-----------------------------------------------------------------------------*/

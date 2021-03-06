#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief Creates system matrix \f$ [A] \f$.
*
*  It is split into separete parts: discretizaton of innertial and diffusive
*  parts, boundary conditions and. 
*  Finally, it also computes the inverse of the system matrix diagonal.
*******************************************************************************/
void EnthalpyFD::create_system(const Scalar * diff_eddy) {
                             
  /*-----------------------------+ 
  |  create system in the core,  |
  |  correct at the boundaries.  |
  +-----------------------------*/
  create_system_innertial();
  create_system_diffusive(diff_eddy);
  create_system_bnd();
#if 0
  std::cout<<"create_system: "<<A.t[1][1][1]<<" "<<A.b[1][1][2]<<"\n";
  exit(0);
#endif

  /*----------------------------------------------+ 
  |  compute invlerse of the central coefficient  |
  +----------------------------------------------*/
  for_ijk(i,j,k) 
    A.ci[i][j][k] = 1.0 / A.c[i][j][k];

  A.ci.exchange();
}

/*-----------------------------------------------------------------------------+
 '$Id: enthalpyfd_create_system.cpp,v 1.4 2015/02/16 08:29:02 sato Exp $'/
+-----------------------------------------------------------------------------*/

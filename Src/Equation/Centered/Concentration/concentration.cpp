#include "concentration.h"

/***************************************************************************//**
*  Initializes parent (Centered), inserts boundary conditions and 
*  discretizes system matrix.
*******************************************************************************/
Concentration::Concentration(const Scalar & PHI, 
                             const Scalar & F,
                             const Vector & U, 
                             Times & T, 
                             Krylov * S,
                             Matter * f) :
/*---------------------+ 
|  initialize parent   |
+---------------------*/
  Centered( PHI.domain(), PHI, F, & U, T, f, NULL, S ) /* NULL is for solid */
{ 
  assert(PHI.domain() == F.domain());
  assert(PHI.domain() == U.domain());

  phi.bnd_update(); // must be called before discretization because of 
                    // dirichlet boundary condition etc.

  discretize();
}	

/******************************************************************************/
Concentration::~Concentration() {
}	

/*-----------------------------------------------------------------------------+
 '$Id: concentration.cpp,v 1.12 2017/10/26 09:22:29 bures_l Exp $'/
+-----------------------------------------------------------------------------*/

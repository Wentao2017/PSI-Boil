#include "vof.h"
#include <cmath>

/******************************************************************************/
VOF::VOF(const Scalar & PHI, 
                   const Scalar & F,
                   const Scalar & K,
                   const real & con, 
                   const real & den, 
                   const Vector & U, 
                   Times & T,
                   Krylov * S) :
/*---------------------+ 
|  initialize parent   | NULL is for solid
+---------------------*/
  jelly( *PHI.domain() ),
  Centered( PHI.domain(), PHI, F , & U, T, &jelly, NULL, S ),
  kappa( &K ),
  clr( *PHI.domain() ),
  nx( *PHI.domain() ),
  ny( *PHI.domain() ),
  nz( *PHI.domain() ),
  nmag( *PHI.domain() ),
  gpx( *PHI.domain() ),
  gpy( *PHI.domain() ),
  gpz( *PHI.domain() ),
  gpxn( *PHI.domain() ),
  gpyn( *PHI.domain() ),
  gpzn( *PHI.domain() ),
  clrn( *PHI.domain() )


/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
  clr     = phi.shape();
  nx      = phi.shape();
  ny      = phi.shape();
  nz      = phi.shape();
  nmag    = phi.shape();
  gpx     = phi.shape();
  gpy     = phi.shape();
  gpz     = phi.shape();
  gpxn    = phi.shape();
  gpyn    = phi.shape();
  gpzn    = phi.shape();
  clrn    = phi.shape();

  assert(PHI.domain() == F.domain());

  pi = acos(-1.0);
  tanfac = 0.90;
  dxmin=std::min(phi.dxc(1),std::min(phi.dyc(1),phi.dzc(1)));
  boil::cart.min_real(&dxmin);
  ww=1.0*dxmin;

  epsnorm=1.0e-12;

  alloc3d(& iflag, phi.ni(), phi.nj(), phi.nk());

  discretize();
//  init();
}	

/******************************************************************************/
VOF::~VOF() {
}	

/*-----------------------------------------------------------------------------+
 '$Id: colorcip.cpp,v 1.3 2009/11/12 12:15:47 sato Exp $'/
+-----------------------------------------------------------------------------*/


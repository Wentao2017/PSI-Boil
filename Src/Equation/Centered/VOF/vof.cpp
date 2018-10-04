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
  clrn( *PHI.domain() ),
  vma( *PHI.domain() ),
  vmb( *PHI.domain() ),
  vmc( *PHI.domain() ),
  alpha( *PHI.domain() ),
  vm1( *PHI.domain() ),
  vm2( *PHI.domain() ),
  vm3( *PHI.domain() ),
  vm12( *PHI.domain() ),
  w( *PHI.domain() ),
  v1( *PHI.domain() ),
  v3( *PHI.domain() ),
  a0( *PHI.domain() ),
  a1( *PHI.domain() ),
  a2( *PHI.domain() ),
  q0( *PHI.domain() ),
  sp( *PHI.domain() ),
  th( *PHI.domain() ),
  absgu( *PHI.domain() ),
  ra( *PHI.domain() ),
  qa( *PHI.domain() ),
  alpha_tmp( *PHI.domain() ),
  vma_tmp( *PHI.domain() ),
  vmb_tmp( *PHI.domain() ),
  vmc_tmp( *PHI.domain() ),
  a( *PHI.domain() ),
  vv( *PHI.domain() ),
  flux_x( *PHI.domain() )

/*------------------------------------------------------+
|  this constructor is called only at the finest level  |
+------------------------------------------------------*/
{ 
  clr       = phi.shape();
  nx        = phi.shape();
  ny        = phi.shape();
  nz        = phi.shape();
  nmag      = phi.shape();
  gpx       = phi.shape();
  gpy       = phi.shape();
  gpz       = phi.shape();
  gpxn      = phi.shape();
  gpyn      = phi.shape();
  gpzn      = phi.shape();
  clrn      = phi.shape();
  vma       = phi.shape();
  vmb       = phi.shape();
  vmc       = phi.shape();
  alpha     = phi.shape();
  vm1       = phi.shape();
  vm2       = phi.shape();
  vm3       = phi.shape();
  vm12      = phi.shape();
  w         = phi.shape();
  v1        = phi.shape();
  v3        = phi.shape(); 
  a0        = phi.shape();
  a1        = phi.shape();
  a2        = phi.shape();
  q0        = phi.shape();
  sp        = phi.shape();
  th        = phi.shape();
  absgu     = phi.shape();
  ra        = phi.shape();
  qa        = phi.shape();
  alpha_tmp = phi.shape();
  vma_tmp   = phi.shape();
  vmb_tmp   = phi.shape();
  vmc_tmp   = phi.shape();
  a         = phi.shape();
  vv        = phi.shape();
  flux_x    = phi.shape();

  assert(PHI.domain() == F.domain());

  pi = acos(-1.0);
  tanfac = 0.90;
  dxmin=std::min(phi.dxc(1),std::min(phi.dyc(1),phi.dzc(1)));
  boil::cart.min_real(&dxmin);
  ww=1.0*dxmin;

  epsnorm=1.0e-12;

/* set initial value */
  nlayer=16;

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


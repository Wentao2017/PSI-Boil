#ifndef VOF_H
#define VOF_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"

#define USE_TAN

///////////
//       //
//  VOF  //
//       //
///////////
class VOF : public Centered {
  public:
    VOF(const Scalar & phi,
        const Scalar & f,
        const Scalar & kappa,
        const real & con, 
        const real & den,
        const Vector & u, 
        Times & t,
        Krylov * S);
    ~VOF();

    void advance();
    void tension(Vector * vec, const Matter matt);
    void totalvol();
    void front_minmax();
    void init();
    void plic_advance();

    // getter for front_minmax
    real get_xminft() { return(xminft);};
    real get_xmaxft() { return(xmaxft);};
    real get_yminft() { return(yminft);};
    real get_ymaxft() { return(ymaxft);};
    real get_zminft() { return(zminft);};
    real get_zmaxft() { return(zmaxft);};

  protected:
    void curv(const Scalar & g);
    void bdcurv(const Scalar & g, const real & v);
    void insert_bc(const Scalar & g);
    void gradphi(const Scalar & g);

    Scalar nx,ny,nz,nmag;/* normal to interface */
    Scalar clr,clrn;     /* color function */
    Scalar gpx,gpy,gpz,gpxn,gpyn,gpzn;
    Scalar kappa;        /* curvature */
    Scalar stmp;
    Matter jelly;   /* virtual fluid for level set transport */
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real pi,tanfac,theta;
    real dxmin,ww;
    real epsnorm;

    int *** iflag;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: colorcip.h,v 1.6 2009/11/12 12:24:07 sato Exp $'/
+-----------------------------------------------------------------------------*/


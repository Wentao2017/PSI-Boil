#ifndef VOF_H
#define VOF_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"

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
        const Vector & u, 
        Times & t,
        Krylov * S);
    ~VOF();

    void new_time_step(){};
    void advance();
    void curvature();
    void tension(Vector * vec, const Matter matt);
    void totalvol();
    void front_minmax();
    void init(){};

    // getter for front_minmax
    real get_xminft() { return(xminft);};
    real get_xmaxft() { return(xmaxft);};
    real get_yminft() { return(yminft);};
    real get_ymaxft() { return(ymaxft);};
    real get_zminft() { return(zminft);};
    real get_zmaxft() { return(zmaxft);};

  protected:
    void advance_x();
    void advance_y();
    void advance_z();
    void bdcurv(const Scalar & g, const real & v);
    void cal_fs();
    void cal_fs2();
    void curv_HF();
    void ext_fs();
    real extract_alpha(const int i, const int j, const int k);
    void insert_bc(const Scalar & g);
    void gradphi(const Scalar & g);
    void gradphic(const Scalar & g);
    void insert_bc_gradphic(const Scalar & g);
    void insert_bc_norm_cc(const Scalar & g);
    void insert_bc_norm();
    void norm_cc(const Scalar & g);
    void normalize(real & r1, real & r2, real & r3);
    real calc_alpha(real & r1, real & r2, real & r3, real & r4);
    real calc_v(real r1, real r2, real r3, real r4);
    void selectMax(const real r1, const real r2, const real r3,
                   const real r4, const real r5, const real r6,
                   const real r7, const real r8, const real r9,
                   const int i1,  const int i2,  const int i3);
    void set_iflag();
    void insert_bc_flag(ScalarInt & g, const bool b);

    void norm_cc_imin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_imax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_jmax(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmin(const Scalar &g, const int i,const int j, const int k);
    void norm_cc_kmax(const Scalar &g, const int i,const int j, const int k);


    Scalar nx,ny,nz,nmag;/* normal to interface */
    Scalar clr,clrn;     /* color function */
    Scalar kappa;        /* curvature */
    Scalar stmp;
    Scalar fsx,fsy,fsz;
    ScalarInt iflag,iflagx,iflagy,iflagz;

    Matter jelly;   /* virtual fluid for level set transport */
    real xminft,xmaxft,yminft,ymaxft,zminft,zmaxft; /* xyz min&max of front */
    real pi,theta;
    real dxmin,ww;
    real epsnorm;
    real phisurf;
    
    int nlayer, n_ext_fs;
    //int *** iflag;
};	
#endif

/*-----------------------------------------------------------------------------+
 '$Id: colorcip.h,v 1.6 2009/11/12 12:24:07 sato Exp $'/
+-----------------------------------------------------------------------------*/


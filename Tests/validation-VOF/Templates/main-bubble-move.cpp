/*----------------------+
|                       |
|  set lcnd to 0.5 ...  |
|                       |
+----------------------*/
#include "Include/psi-boil.h"

/* boundary conditions */
const int NX= 32;
const int NY= NX*2;
const int NZ= NX;

const real RB = 0.0025;           
const real LX = RB*4.0;           
const real LY = LX*NY/NX;
const real LZ = LX; 

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*----------+
  |  grid(s)  |
  +----------*/

//  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), NX, Periodic::yes());
//  Grid1D gy( Range<real>(-LY/2.0, LY/2.0), NY, Periodic::yes());
//  Grid1D gz( Range<real>(-LZ/2.0, LZ/2.0), NZ, Periodic::yes());

  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), Range<real>(LZ/(2.0 * NZ),LZ/(2.0 * NZ)), NX, Periodic::yes());
  Grid1D gy( Range<real>(-LY/2.0, LY/2.0), Range<real>(LZ/(2.0 * NZ),LZ/(2.0 * NZ)), NY, Periodic::yes());
  Grid1D gz( Range<real>(-LZ/2.0, LZ/2.0), Range<real>(LZ/(2.0 * NZ),LZ/(2.0 * NZ)), NZ, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*-------------------+
  |  time-integration  |
  +-------------------*/

  const int  ndt = 60;
  const real dt  = 0.25 * LX / real(NX);
  const int  nint = 1;

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d);
  Scalar c  (d), g  (d); // concentration
  Scalar kappa(d);

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  c.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
 // c.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 0.0 ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );

#if 0
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  }
#endif



  /*--------------------+
  |  initial condition  |
  +--------------------*/

  Comp m=Comp::u();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=-1.0;

  m=Comp::v();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=-1.0;

  m=Comp::w();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=-1.0;

  uvw.exchange_all();


  for_vijk(c,i,j,k)
      c[i][j][k] = 1.0;

  const real radius=RB;
  const real zcent =0;
  for_vijk(c,i,j,k) {
    //if (c.zc(k)>0.0) {
    real dist=pow(c.xc(i),2.0)+pow(c.yc(j),2.0)+pow((c.zc(k)-zcent),2.0);
    if (dist<pow(radius*0.75,2.0)) {
      c[i][j][k]=0.0;
    } else if(dist<pow(radius*1.25,2.0)) {
      int mm=8;
      real x0=d.xn(i);
      real y0=d.yn(j);
      real z0=d.zn(k);
      real ddx=d.dxc(i)/real(mm);
      real ddy=d.dyc(j)/real(mm);
      real ddz=d.dzc(k)/real(mm);
      int itmp=0;
      for (int ii=0; ii<mm; ii++){
        for (int jj=0; jj<mm; jj++){
          for (int kk=0; kk<mm; kk++){
            real xxc=x0+0.5*ddx+real(ii)*ddx;
            real yyc=y0+0.5*ddy+real(jj)*ddy;
            real zzc=z0+0.5*ddz+real(kk)*ddz;
            real dist=pow(xxc,2.0)+pow(yyc,2.0)+pow(zzc-zcent,2.0);
            if (dist>pow(radius,2.0)){
                itmp=itmp+1;
            }
          }
        }
      }
      c[i][j][k]=real(itmp)/real(mm*mm*mm);
    //}
    }
  }
  c.exchange_all();
  c.bnd_update();
  boil::plot->plot(uvw,c, "uvw-c", 0); 
      
  Times time(ndt, dt); /* ndt, dt */

  VOF conc(c, g, kappa, 1.0, 1.0, uvw, time, solver);	

  for(time.start(); time.end(); time.increase()) {

  boil::oout << "########################" << boil::endl;
  boil::oout << "#                       " << boil::endl;
  boil::oout << "# TIME:      " << time.current_time() << boil::endl;
  boil::oout << "#                       " << boil::endl;
  boil::oout << "# TIME STEP: " << time.current_step()
             << "/"             << time.total_steps() << boil::endl;
  boil::oout << "#                       " << boil::endl;
  boil::oout << "########################" << boil::endl;

  
  conc.advance();
  conc.totalvol();

  if(time.current_step() % nint == 0) {
      boil::plot->plot(uvw, c, "uvw-c",  time.current_step());
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
  
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-zeppelin.cpp,v 1.2 2012/09/13 08:13:33 niceno Exp $'/
+-----------------------------------------------------------------------------*/

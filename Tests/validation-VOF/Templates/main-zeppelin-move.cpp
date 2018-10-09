/*----------------------+
|                       |
|  set lcnd to 0.5 ...  |
|                       |
+----------------------*/
#include "Include/psi-boil.h"
#include <fstream>

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
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-LX/2.0, LX/2.0), NX, Periodic::no());
  Grid1D gy( Range<real>(-LY/2.0, LY/2.0), NY, Periodic::no());
  Grid1D gz( Range<real>(-LZ/2.0, LZ/2.0), NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 20;
  //const int ndt = 1;
  const int nint = 1;
  const real dt  = 0.25 * LX / real(NX);
  Times time(ndt, dt);
  
  OPR(  NX );
  OPR(  dt );
  OPR( ndt );
  
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
  c.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 0.0 ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  /*--------------------+
  |  initial condition  |
  +--------------------*/
  Comp m=Comp::u();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=1.0;

  m=Comp::v();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=0.0;

  m=Comp::w();
  for_avmijk(uvw,m,i,j,k)
    uvw[m][i][j][k]=0.0;

  uvw.exchange_all();

  /*---------+
  |  circle  |
  +---------*/
  for_vijk(c,i,j,k) {
    real dist = sqrt( c.xc(i)*c.xc(i) + 0.25*c.yc(j)*c.yc(j) + c.zc(k)*c.zc(k) ); 
      if(dist < RB)
        c[i][j][k] = 1.0; /* gas   */
      else
        c[i][j][k] = 0.0; /* liquid*/
  }

  c.exchange();

  ColorFunction conc  (c, g, 1.0, 1.0, uvw, time, solver);

  for(time.start(); time.end(); time.increase()) {

  boil::oout << "########################" << boil::endl;
  boil::oout << "#                       " << boil::endl;
  boil::oout << "# TIME:      " << time.current_time() << boil::endl;
  boil::oout << "#                       " << boil::endl;
  boil::oout << "# TIME STEP: " << time.current_step()
             << "/"             << time.total_steps() << boil::endl;
  boil::oout << "#                       " << boil::endl;
  boil::oout << "########################" << boil::endl;

    /*---------------------------+
    |  fully explicit with conc  |
    +---------------------------*/
    conc.new_time_step();
    //conc.convection();
    conc.advance();
    conc.sharpen();
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

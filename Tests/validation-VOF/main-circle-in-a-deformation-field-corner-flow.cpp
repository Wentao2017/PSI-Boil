#include "Include/psi-boil.h"
//#define ColorFunc

/* domain dimensions (given by problem) */
const real LX =   1.0;
const real LZ =   0.05; //0.004;

/* computed parameters */
const int NX = 100;
const int NZ = 4;

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
  Grid1D gx( Range<real>( -0.5*LX,0.5*LX), NX, Periodic::yes() );
  Grid1D gz( Range<real>( -0.5*LZ,0.5*LZ), NZ, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gz);

  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d), xyz(d); // vel
  Scalar p  (d), f  (d); // p.
  Scalar c  (d), g  (d), kappa(d); // concentration

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  }
  
  p.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );

  /* copy b.c. from p */
  
  c.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );

  g = c.shape();


  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt =  200;
  const int  nint =  50;
  const real dt  = 0.0025;

  Times time(ndt, dt);

  OPR(  NX );
  OPR(  dt );
  OPR( ndt );

  for_vijk(c,i,j,k) 
    c[i][j][k] = 0.0;

  const real radius=0.15;
  const real xcent=0.0;
  const real ycent=0.0;
  for_vijk(c,i,j,k) {
    real dist=pow(c.xc(i)-xcent,2.0)+pow(c.yc(j)-ycent,2.0);
    if (dist<pow(radius*0.75,2.0)) {
      c[i][j][k]=1.0;
    } else if(dist<pow(radius*1.25,2.0)) {
      int mm=10;
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
            real dist=pow(xxc-xcent,2.0)+pow(yyc-ycent,2.0);
            if (dist<pow(radius,2.0)){
              itmp=itmp+1;
            }
          }
        }
      }
      c[i][j][k]=real(itmp)/real(mm*mm*mm);
    }
  }
  c.exchange_all();

  VOF conc  (c, g, kappa, 1.0, 1.0, uvw, time, solver);
  boil::plot->plot(c,"c", 0);


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
#if 1
    conc.new_time_step();

    real pi=acos(-1.0);
    real period=8.0;
    real t=time.current_time();
    Comp m = Comp::u();
    for_avmijk(uvw,m,i,j,k){
      real x=uvw.xc(m,i);
      real y=uvw.yc(m,j);
      uvw[m][i][j][k] =  2.0*x;
    }
    m = Comp::v();
    for_avmijk(uvw,m,i,j,k){
      real x=uvw.xc(m,i);
      real y=uvw.yc(m,j);
      uvw[m][i][j][k] = -2.0*y;
    }

    uvw.exchange_all();
#endif 
   
    conc.advance();
    conc.totalvol();

    if(time.current_step() % (nint)==0) {
      boil::plot->plot(uvw, c, "uvw-c",  time.current_step());
    }
    if(time.current_step()==1) {
      boil::plot->plot(uvw, c, "uvw-c", 0);
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}	
/*-----------------------------------------------------------------------------+
 '$Id: main-corner-flow.cpp,v 1.4 2018/09/26 10:09:00 sato Exp $'/
+-----------------------------------------------------------------------------*/

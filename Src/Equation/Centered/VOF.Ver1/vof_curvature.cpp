#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::curvature() {
/***************************************************************************//**
*  \brief Calculate curvature using height function.
*     output: kappa
*******************************************************************************/

  /*--------------------------+
  |  flaggin interface cells  |
  +--------------------------*/
  set_iflag();

  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
  norm_cc(phi);

  /*----------------------------------+
  |  calculate free surface position  |
  +----------------------------------*/
  cal_fs();
  //cal_fs2();

  /*------------------------------------+
  |  extrapolate free surface position  |
  +------------------------------------*/
  ext_fs();

  /* interfacial cells: iflag=1 */
#if 0
  for_ijk(i,j,k) {
    real sum_phi=0.0;
    for(int ii=-1; ii<=1; ii++) {
      for(int jj=-1; jj<=1; jj++) {
        for(int kk=-1; kk<=1; kk++) {
          sum_phi += phi[i+ii][j+jj][k+kk];
        }
      }
    }
    if ((boil::micro<sum_phi)&&(sum_phi<27.0-boil::micro)) {
      iflag[i][j][k]=1;
    }
  }
  iflag.exchange();
  // crude code: need boundary treatment for iflag //
#endif

  /* curvature calculation */
  kappa=0.0;
  for_ijk(i,j,k) {
    //if (iflag[i][j][k]==1) {
#if 0
      if(i==27&&j==23&&k==23) {
        cout<<"iflag= "<<iflag[i][j][k]<<"\n";
      }
#endif

    if (abs(iflag[i][j][k])<=2) {

#if 0
      if(i==27&&j==23&&k==23) {
        cout<<"iflag= "<<iflag[i][j][k]<<"\n";
      }
#endif

      real dfsx = fsx[i][j][k]-phi.xc(i);
      real dfsy = fsy[i][j][k]-phi.yc(j);
      real dfsz = fsz[i][j][k]-phi.zc(k);

      int n_in=0;
      if (fabs(dfsx)<0.5*phi.dxc(i)) n_in++;  //fsx is in the cell
      if (fabs(dfsy)<0.5*phi.dyc(j)) n_in++;  //fsy is in the cell
      if (fabs(dfsz)<0.5*phi.dzc(k)) n_in++;  //fsz is in the cell

      //if (i==13&&j==1&&k==50) {
      //  std::cout<<dfsx<<" "<<dfsy<<" "<<dfsz<<" "<<n_in<<"\n";
      //  exit(0);
      //}
      int dirMax=0;
      if (n_in>=2) {

        //std::cout<<"n_in>2: "<<n_in<<" "<<i<<" "<<j<<" "<<k<<"\n";
        real nxx = (phi[i+1][j][k]-phi[i-1][j][k])/(dxw(i)+dxe(i));
        real nyy = (phi[i][j+1][k]-phi[i][j-1][k])/(dys(j)+dyn(j));
        real nzz = (phi[i][j][k+1]-phi[i][j][k-1])/(dzb(k)+dzt(k));

        if (fabs(nxx)<fabs(nyy)) {
          if (fabs(nyy)<fabs(nzz)) {
            dirMax=3;
          } else {
            dirMax=2;
          }
        } else {
          if (fabs(nxx)<fabs(nzz)) {
            dirMax=3;
          } else {
            dirMax=1;
          }
        }
#if 0
        if (i==34&&j==68&&k==57) {
          cout<<"dir= "<<dirMax<<" "<<n_in<<" "<<nxx<<" "<<nyy<<" "<<nzz<<"\n";
        }
#endif

      } else {

        if (fabs(dfsx)>fabs(dfsy)) {
          if (fabs(dfsy)>fabs(dfsz)) {
            dirMax=3;
          } else {
            dirMax=2;
          }
        } else {
          if (fabs(dfsx)>fabs(dfsz)) {
            dirMax=3;
          } else {
            dirMax=1;
          }
        }
      }

#if 0
      if (i==27&&j==1&&k==57) {
        cout<<"dir= "<<dirMax<<" "<<n_in<<"\n";
      }
#endif
 
      if (dirMax==1) {
        real hmm, hcm, hpm, hmc, hcc, hpc, hmp, hcp, hpp;
#if 0
        hmm = phi[i-1][j-1][k-1]*phi.dxc(i-1)
            + phi[i  ][j-1][k-1]*phi.dxc(i)
            + phi[i+1][j-1][k-1]*phi.dxc(i+1);
        hcm = phi[i-1][j  ][k-1]*phi.dxc(i-1)
            + phi[i  ][j  ][k-1]*phi.dxc(i)
            + phi[i+1][j  ][k-1]*phi.dxc(i+1);
        hpm = phi[i-1][j+1][k-1]*phi.dxc(i-1)
            + phi[i  ][j+1][k-1]*phi.dxc(i)
            + phi[i+1][j+1][k-1]*phi.dxc(i+1);
        hmc = phi[i-1][j-1][k  ]*phi.dxc(i-1)
            + phi[i  ][j-1][k  ]*phi.dxc(i)
            + phi[i+1][j-1][k  ]*phi.dxc(i+1);
        hcc = phi[i-1][j  ][k  ]*phi.dxc(i-1)
            + phi[i  ][j  ][k  ]*phi.dxc(i)
            + phi[i+1][j  ][k  ]*phi.dxc(i+1);
        hpc = phi[i-1][j+1][k  ]*phi.dxc(i-1)
            + phi[i  ][j+1][k  ]*phi.dxc(i)
            + phi[i+1][j+1][k  ]*phi.dxc(i+1);
        hmp = phi[i-1][j-1][k+1]*phi.dxc(i-1)
            + phi[i  ][j-1][k+1]*phi.dxc(i)
            + phi[i+1][j-1][k+1]*phi.dxc(i+1);
        hcp = phi[i-1][j  ][k+1]*phi.dxc(i-1)
            + phi[i  ][j  ][k+1]*phi.dxc(i)
            + phi[i+1][j  ][k+1]*phi.dxc(i+1);
        hpp = phi[i-1][j+1][k+1]*phi.dxc(i-1)
            + phi[i  ][j+1][k+1]*phi.dxc(i)
            + phi[i+1][j+1][k+1]*phi.dxc(i+1);
#else
        real xcc = phi.xc(i);
        real hmax = phi.dxc(i)*real(n_ext_fs+1);  // crude code
        hmm = copysign(1.0,phi[i][j-1][k-1]-0.5)
             * min(fabs(fsx[i][j-1][k-1]-xcc),hmax);
        hcm = copysign(1.0,phi[i][j  ][k-1]-0.5)
             * min(fabs(fsx[i][j  ][k-1]-xcc),hmax);
        hpm = copysign(1.0,phi[i][j+1][k-1]-0.5)
             * min(fabs(fsx[i][j+1][k-1]-xcc),hmax);
        hmc = copysign(1.0,phi[i][j-1][k  ]-0.5)
             * min(fabs(fsx[i][j-1][k  ]-xcc),hmax);
        hcc = copysign(1.0,phi[i][j  ][k  ]-0.5)
             * min(fabs(fsx[i][j  ][k  ]-xcc),hmax);
        hpc = copysign(1.0,phi[i][j+1][k  ]-0.5)
             * min(fabs(fsx[i][j+1][k  ]-xcc),hmax);
        hmp = copysign(1.0,phi[i][j-1][k+1]-0.5)
             * min(fabs(fsx[i][j-1][k+1]-xcc),hmax);
        hcp = copysign(1.0,phi[i][j  ][k+1]-0.5)
             * min(fabs(fsx[i][j  ][k+1]-xcc),hmax);
        hpp = copysign(1.0,phi[i][j+1][k+1]-0.5)
             * min(fabs(fsx[i][j+1][k+1]-xcc),hmax);
        if(hmm+hcm+hpm+hmc+hcc+hpc+hmp+hcp+hpp>boil::zetta) {
          std::cout<<"Error:"<<i<<" "<<j<<" "<<k<<" "<<phi[i][j][k]<<" "
                   <<"dir= "<<dirMax<<" "<<dfsx<<" "<<dfsy<<" "<<dfsz<<"\n";
          std::cout<<hmm<<" "<<hcm<<" "<<hpm<<"\n";
          std::cout<<hmc<<" "<<hcc<<" "<<hpc<<"\n";
          std::cout<<hmp<<" "<<hcp<<" "<<hpp<<"\n";
          exit(0);
        }
#if 0
        if (i==12&&j==1&&k==48) {
          std::cout<<"dir1:"<<i<<" "<<j<<" "<<k<<" "<<phi[i][j][k]<<" "
          //         <<nxx<<" "<<nyy<<" "<<nzz<<"\n";
                   <<dirMax<<" "<<dfsx<<" "<<dfsy<<" "<<dfsz<<"\n";
          std::cout<<hmm<<" "<<hcm<<" "<<hpm<<"\n";
          std::cout<<hmc<<" "<<hcc<<" "<<hpc<<"\n";
          std::cout<<hmp<<" "<<hcp<<" "<<hpp<<"\n";
          exit(0);
        }
#endif
#endif
        real hy  = (hpc-hmc)/(dys(j)+dyn(j));
        real hz  = (hcp-hcm)/(dzb(k)+dzt(k));
        real hyy = (hpc-2.0*hcc+hmc)/(phi.dyc(j)*phi.dyc(j));
        real hzz = (hcp-2.0*hcc+hcm)/(phi.dzc(k)*phi.dzc(k));
        real hyz = (hpp-hpm-hmp+hmm)
                 / (4.0*phi.dyc(j)*phi.dzc(k));


        real curv = -1.0
                    * (hyy + hzz + hyy*hz*hz + hzz*hy*hy + 2.0*hyz*hy*hz)
                    / pow(1.0 + hy*hy + hz*hz, 1.5);
        real curv_max = 1.0/phi.dyc(j) + 1.0/phi.dzc(k);
        real curv_min = -curv_max;

        kappa[i][j][k] = max(curv_min,min(curv_max,curv));
#if 0
        if (i==27&&j==1&&k==57) {
          cout<<"hmax= "<<hmax<<" "<<kappa[i][j][k]<<"\n";
          cout<<"hx= "<<hy<<" "<<hz<<" "<<hyy<<" "<<hzz<<" "<<hyz<<"\n";
          cout<<"hmm= "<<hmm<<" "<<hcm<<" "<<hpm<<"\n";
          cout<<"hmc= "<<hmc<<" "<<hcc<<" "<<hpc<<"\n";
          cout<<"hmp= "<<hmp<<" "<<hcp<<" "<<hpp<<"\n";
          cout<<"fsx= "<<fsx[i][j][k]<<" "<<fsx[i][j][k-1]<<"\n";
        }
#endif

      } else if (dirMax==2) {
        real hmm, hcm, hpm, hmc, hcc, hpc, hmp, hcp, hpp;
        real ycc = phi.yc(j);
        real hmax = phi.dyc(j)*real(n_ext_fs+1);  // crude code
        hmm = copysign(1.0,phi[i-1][j][k-1]-0.5)
             * min(fabs(fsy[i-1][j][k-1]-ycc),hmax);
        hcm = copysign(1.0,phi[i  ][j][k-1]-0.5)
             * min(fabs(fsy[i  ][j][k-1]-ycc),hmax);
        hpm = copysign(1.0,phi[i+1][j][k-1]-0.5)
             * min(fabs(fsy[i+1][j][k-1]-ycc),hmax);
        hmc = copysign(1.0,phi[i-1][j][k  ]-0.5)
             * min(fabs(fsy[i-1][j][k  ]-ycc),hmax);
        hcc = copysign(1.0,phi[i  ][j][k  ]-0.5)
             * min(fabs(fsy[i  ][j][k  ]-ycc),hmax);
        hpc = copysign(1.0,phi[i+1][j][k  ]-0.5)
             * min(fabs(fsy[i+1][j][k  ]-ycc),hmax);
        hmp = copysign(1.0,phi[i-1][j][k+1]-0.5)
             * min(fabs(fsy[i-1][j][k+1]-ycc),hmax);
        hcp = copysign(1.0,phi[i  ][j][k+1]-0.5)
             * min(fabs(fsy[i  ][j][k+1]-ycc),hmax);
        hpp = copysign(1.0,phi[i+1][j][k+1]-0.5)
             * min(fabs(fsy[i+1][j][k+1]-ycc),hmax);
        if(hmm+hcm+hpm+hmc+hcc+hpc+hmp+hcp+hpp>boil::zetta) {
          std::cout<<"Error:"<<i<<" "<<j<<" "<<k<<" "<<phi[i][j][k]<<" "
          //         <<nxx<<" "<<nyy<<" "<<nzz<<"\n";
                   <<"dir= "<<dirMax<<" "<<dfsx<<" "<<dfsy<<" "<<dfsz<<"\n";
          std::cout<<hmm<<" "<<hcm<<" "<<hpm<<"\n";
          std::cout<<hmc<<" "<<hcc<<" "<<hpc<<"\n";
          std::cout<<hmp<<" "<<hcp<<" "<<hpp<<"\n";
          exit(0);
        }
        real hx  = (hpc-hmc)/(dxw(i)+dxe(i));
        real hz  = (hcp-hcm)/(dzb(k)+dzt(k));
        real hxx = (hpc-2.0*hcc+hmc)/(phi.dxc(i)*phi.dxc(i));
        real hzz = (hcp-2.0*hcc+hcm)/(phi.dzc(k)*phi.dzc(k));
        real hxz = (hpp-hpm-hmp+hmm)
                 / (4.0*phi.dxc(i)*phi.dzc(k));

        real curv = -1.0
                  * (hxx + hzz + hxx*hz*hz + hzz*hx*hx + 2.0*hxz*hx*hz)
                  / pow(1.0 + hx*hx + hz*hz, 1.5);
        real curv_max = 1.0/phi.dxc(i) + 1.0/phi.dzc(k);
        real curv_min = -curv_max;

        kappa[i][j][k] = max(curv_min,min(curv_max,curv));

#if 0
        if (i==27&&j==23&&k==23) {
          cout<<"hmax= "<<hmax<<" "<<kappa[i][j][k]<<"\n";
          cout<<"hx= "<<hx<<" "<<hz<<" "<<hxx<<" "<<hzz<<" "<<hxz<<"\n";
          cout<<"hmm= "<<hmm<<" "<<hcm<<" "<<hpm<<"\n";
          cout<<"hmc= "<<hmc<<" "<<hcc<<" "<<hpc<<"\n";
          cout<<"hmp= "<<hmp<<" "<<hcp<<" "<<hpp<<"\n";
          cout<<"i-1,k-1= "<<phi[i-1][j][k-1]<<" "<<fsy[i-1][j][k-1]-ycc<<"\n";
        }
#endif

      } else if (dirMax==3) {
        // calculate height
        real hmm, hcm, hpm, hmc, hcc, hpc, hmp, hcp, hpp;
#if 0
        hmm = phi[i-1][j-1][k-1]*phi.dzc(k-1)
            + phi[i-1][j-1][k  ]*phi.dzc(k)
            + phi[i-1][j-1][k+1]*phi.dzc(k+1);
        hcm = phi[i  ][j-1][k-1]*phi.dzc(k-1)
            + phi[i  ][j-1][k  ]*phi.dzc(k)
            + phi[i  ][j-1][k+1]*phi.dzc(k+1);
        hpm = phi[i+1][j-1][k-1]*phi.dzc(k-1)
            + phi[i+1][j-1][k  ]*phi.dzc(k)
            + phi[i+1][j-1][k+1]*phi.dzc(k+1);
        hmc = phi[i-1][j  ][k-1]*phi.dzc(k-1)
            + phi[i-1][j  ][k  ]*phi.dzc(k)
            + phi[i-1][j  ][k+1]*phi.dzc(k+1);
        hcc = phi[i  ][j  ][k-1]*phi.dzc(k-1)
            + phi[i  ][j  ][k  ]*phi.dzc(k)
            + phi[i  ][j  ][k+1]*phi.dzc(k+1);
        hpc = phi[i+1][j  ][k-1]*phi.dzc(k-1)
            + phi[i+1][j  ][k  ]*phi.dzc(k)
            + phi[i+1][j  ][k+1]*phi.dzc(k+1);
        hmp = phi[i-1][j+1][k-1]*phi.dzc(k-1)
            + phi[i-1][j+1][k  ]*phi.dzc(k)
            + phi[i-1][j+1][k+1]*phi.dzc(k+1);
        hcp = phi[i  ][j+1][k-1]*phi.dzc(k-1)
            + phi[i  ][j+1][k  ]*phi.dzc(k)
            + phi[i  ][j+1][k+1]*phi.dzc(k+1);
        hpp = phi[i+1][j+1][k-1]*phi.dzc(k-1)
            + phi[i+1][j+1][k  ]*phi.dzc(k)
            + phi[i+1][j+1][k+1]*phi.dzc(k+1);
#else
        real zcc = phi.zc(k);
        real hmax = phi.dzc(k)*real(n_ext_fs+1);
        hmm = copysign(1.0,phi[i-1][j-1][k]-0.5)
             * min (fabs(fsz[i-1][j-1][k]-zcc),hmax);
        hcm = copysign(1.0,phi[i  ][j-1][k]-0.5)
             * min (fabs(fsz[i  ][j-1][k]-zcc),hmax);
        hpm = copysign(1.0,phi[i+1][j-1][k]-0.5)
             * min (fabs(fsz[i+1][j-1][k]-zcc),hmax);
        hmc = copysign(1.0,phi[i-1][j  ][k]-0.5)
             * min (fabs(fsz[i-1][j  ][k]-zcc),hmax);
        hcc = copysign(1.0,phi[i  ][j  ][k]-0.5)
             * min (fabs(fsz[i  ][j  ][k]-zcc),hmax);
        hpc = copysign(1.0,phi[i+1][j  ][k]-0.5)
             * min (fabs(fsz[i+1][j  ][k]-zcc),hmax);
        hmp = copysign(1.0,phi[i-1][j+1][k]-0.5)
             * min (fabs(fsz[i-1][j+1][k]-zcc),hmax);
        hcp = copysign(1.0,phi[i  ][j+1][k]-0.5)
             * min (fabs(fsz[i  ][j+1][k]-zcc),hmax);
        hpp = copysign(1.0,phi[i+1][j+1][k]-0.5)
             * min (fabs(fsz[i+1][j+1][k]-zcc),hmax);
        //if (fabs(fsz[i-1][j-1][k]-zcc)>hmax) { 
        //  cout<<"hmm: "<<fabs(fsz[i-1][j-1][k]-zcc)<<" "<<hmax<<"\n";
        //  exit(0);
        //}
        //if (fabs(fsz[i+1][j+1][k]-zcc)>hmax) { 
        //  cout<<"hpp: "<<fabs(fsz[i+1][j+1][k]-zcc)<<" "<<hmax<<"\n";
        //  exit(0);
        //}
        if(hmm+hcm+hpm+hmc+hcc+hpc+hmp+hcp+hpp>boil::zetta) {
          std::cout<<"Error:"<<i<<" "<<j<<" "<<k<<" "<<phi[i][j][k]<<" "
          //         <<nxx<<" "<<nyy<<" "<<nzz<<"\n";
                   <<"dir= "<<dirMax<<" "<<dfsx<<" "<<dfsy<<" "<<dfsz<<"\n";
          std::cout<<hmm<<" "<<hcm<<" "<<hpm<<"\n";
          std::cout<<hmc<<" "<<hcc<<" "<<hpc<<"\n";
          std::cout<<hmp<<" "<<hcp<<" "<<hpp<<"\n";
          exit(0);
        }
#endif
        real hx  = (hpc-hmc)/(dxw(i)+dxe(i));
        real hy  = (hcp-hcm)/(dys(j)+dyn(j));
        real hxx = (hpc-2.0*hcc+hmc)/(phi.dxc(i)*phi.dxc(i));
        real hyy = (hcp-2.0*hcc+hcm)/(phi.dyc(j)*phi.dyc(j));
        real hxy = (hpp-hpm-hmp+hmm)
                 / (4.0*phi.dxc(i)*phi.dyc(j));

        real curv = -1.0
                  * (hxx + hyy + hxx*hy*hy + hyy*hx*hx + 2.0*hxy*hx*hy)
                  / pow(1.0 + hx*hx + hy*hy, 1.5);
        real curv_max = 1.0/phi.dxc(i) + 1.0/phi.dyc(j);
        real curv_min = -curv_max;

        kappa[i][j][k] = max(curv_min,min(curv_max,curv));
      }
#if 0
      if(i==27&&j==23&&k==23) {
        cout<<"kappa= "<<kappa[i][j][k]<<" "<<dirMax<<"\n";
      }
#endif
    }
  }
  kappa.exchange();

  // curvature at boundary 
  //bd_curve();  need to implement!!!!


#if 0
  for_aijk(i,j,k) {
    stmp[i][j][k]=real(iflag[i][j][k]);
  }
  boil::plot->plot(phi,kappa,stmp, "phi-kappa-iflag", time->current_step());
  exit(0);
#endif

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: colorcip_curv.cpp,v 1.2 2009/11/12 12:15:48 sato Exp $'/
+-----------------------------------------------------------------------------*/


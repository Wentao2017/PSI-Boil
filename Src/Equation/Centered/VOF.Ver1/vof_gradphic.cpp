#include "vof.h"

/******************************************************************************/
void VOF::gradphic(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate grad(csa)/|grad(csa)| at cell center.
*         Resluts: nx, ny, nz
*******************************************************************************/

  /* cell centered base */
  for_ijk(i,j,k) {
#if 0
      nx[i][j][k] = (sca[i+1][j][k]-sca[i-1][j][k])/(dxw(i)+dxe(i));
      ny[i][j][k] = (sca[i][j+1][k]-sca[i][j-1][k])/(dys(j)+dyn(j));
      nz[i][j][k] = (sca[i][j][k+1]-sca[i][j][k-1])/(dzb(k)+dzt(k));
#endif
#if 0
      nx[i][j][k] = (sca[i+1][j][k]-sca[i-1][j][k]);
      ny[i][j][k] = (sca[i][j+1][k]-sca[i][j-1][k]);
      nz[i][j][k] = (sca[i][j][k+1]-sca[i][j][k-1]);
#endif
#if 1
      real q000, q001, q010, q011, q100, q101, q110, q111;
      q000 = 0.125 * (sca[i-1][j-1][k-1] + sca[i][j-1][k-1]
                    + sca[i-1][j  ][k-1] + sca[i][j  ][k-1]
                    + sca[i-1][j-1][k  ] + sca[i][j-1][k  ]
                    + sca[i-1][j  ][k  ] + sca[i][j  ][k  ]);
      q001 = 0.125 * (sca[i-1][j-1][k  ] + sca[i][j-1][k  ]
                    + sca[i-1][j  ][k  ] + sca[i][j  ][k  ]
                    + sca[i-1][j-1][k+1] + sca[i][j-1][k+1]
                    + sca[i-1][j  ][k+1] + sca[i][j  ][k+1]);
      q100 = 0.125 * (sca[i  ][j-1][k-1] + sca[i+1][j-1][k-1]
                    + sca[i  ][j  ][k-1] + sca[i+1][j  ][k-1]
                    + sca[i  ][j-1][k  ] + sca[i+1][j-1][k  ]
                    + sca[i  ][j  ][k  ] + sca[i+1][j  ][k  ]);
      q101 = 0.125 * (sca[i  ][j-1][k  ] + sca[i+1][j-1][k  ]
                    + sca[i  ][j  ][k  ] + sca[i+1][j  ][k  ]
                    + sca[i  ][j-1][k+1] + sca[i+1][j-1][k+1]
                    + sca[i  ][j  ][k+1] + sca[i+1][j  ][k+1]);
      q010 = 0.125 * (sca[i-1][j  ][k-1] + sca[i][j  ][k-1]
                    + sca[i-1][j+1][k-1] + sca[i][j+1][k-1]
                    + sca[i-1][j  ][k  ] + sca[i][j  ][k  ]
                    + sca[i-1][j+1][k  ] + sca[i][j+1][k  ]);
      q011 = 0.125 * (sca[i-1][j  ][k  ] + sca[i][j  ][k  ]
                    + sca[i-1][j+1][k  ] + sca[i][j+1][k  ]
                    + sca[i-1][j  ][k+1] + sca[i][j  ][k+1]
                    + sca[i-1][j+1][k+1] + sca[i][j+1][k+1]);
      q110 = 0.125 * (sca[i  ][j  ][k-1] + sca[i+1][j  ][k-1]
                    + sca[i  ][j+1][k-1] + sca[i+1][j+1][k-1]
                    + sca[i  ][j  ][k  ] + sca[i+1][j  ][k  ]
                    + sca[i  ][j+1][k  ] + sca[i+1][j+1][k  ]);
      q111 = 0.125 * (sca[i  ][j  ][k  ] + sca[i+1][j  ][k  ]
                    + sca[i  ][j+1][k  ] + sca[i+1][j+1][k  ]
                    + sca[i  ][j  ][k+1] + sca[i+1][j  ][k+1]
                    + sca[i  ][j+1][k+1] + sca[i+1][j+1][k+1]);
#if 0
      nx[i][j][k] = 0.25 * ( q100 - q000 + q110 - q010
                           + q101 - q001 + q111 - q011)/phi.dxc(i);
      ny[i][j][k] = 0.25 * ( q010 - q000 + q110 - q100
                           + q011 - q001 + q111 - q101)/phi.dyc(j);
      nz[i][j][k] = 0.25 * ( q001 - q000 + q101 - q100
                           + q011 - q010 + q111 - q110)/phi.dzc(k);
#else
      nx[i][j][k] = 0.25 * ( q100 - q000 + q110 - q010
                           + q101 - q001 + q111 - q011);
      ny[i][j][k] = 0.25 * ( q010 - q000 + q110 - q100
                           + q011 - q001 + q111 - q101);
      nz[i][j][k] = 0.25 * ( q001 - q000 + q101 - q100
                           + q011 - q010 + q111 - q110);

#endif
#endif
  }

//  std::cout<<"nx75 "<<nx[75][1][1]<<" "<<"nx76 "<<nx[76][1][1]<<"\n";
//  std::cout<<"ny75 "<<ny[75][1][1]<<" "<<"ny76 "<<ny[76][1][1]<<"\n";
//  std::cout<<"nz75 "<<nz[75][1][1]<<" "<<"nz76 "<<nz[76][1][1]<<"\n";

  /* normal vector at adjacent cells next to wall, symmetric and IB */
  insert_bc_gradphic(sca); 

  /* normal vector on boundary plane */
  insert_bc_norm();

  /* normalize */
  for_avijk(sca,i,j,k) {
    normalize(nx[i][j][k],ny[i][j][k],nz[i][j][k]);
  }

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_gradphic.cpp,v 1.7 2015/05/05 15:11:07 sato Exp $'/
+-----------------------------------------------------------------------------*/

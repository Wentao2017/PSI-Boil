#include "colorcip.h"

/******************************************************************************/
void ColorCIP::discretize() {

  /* correct on the boundaries */
  create_system_bnd();   

}

/*-----------------------------------------------------------------------------+
 '$Id: colorcip_discretize.cpp,v 1.2 2009/11/05 09:34:41 sato Exp $'/
+-----------------------------------------------------------------------------*/

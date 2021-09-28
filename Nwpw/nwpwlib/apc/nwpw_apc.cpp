/* nwpw_apc.cpp -
   Author - Eric Bylaska
*/


#include        <iostream>
#include        <cstring>
#include        <cmath>

#include        "nwpw_timing.hpp"
#include        "gdevice.hpp"


#include        "blas.h"

#include        "nwpw_apc.hpp"

namespace pwdft {



/* Constructors */

/*******************************************
 *                                         *
 *            nwpw_apc::nwpw_apc           *
 *                                         *
 *******************************************/
nwpw_apc::nwpw_apc(Ion *myionin, Pneb *mypnebin, Strfac *mystrfacin, Control2& control)
{

   myion    = myionin;
   mypneb   = mypneb;
   mystrfac = mystrfacin;
   apc_on = control.APC_on();

   nga = control.APC_nga();
   ngs = nga*myion->nion;


}


}

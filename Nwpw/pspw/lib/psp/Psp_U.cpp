/* Psp_U.cpp - 
   Author - Eric Bylaska
*/

/**
 * @class Psp_U
 * @brief Class for managing uterm data and calculations.
 *
 * The Psp_U class encapsulates operations and data related to
 * uterm used in electronic structure calculations. It provides
 * methods for handling non-local and local pseudopotentials, semicore
 * corrections, and other properties used in electronic structure calculations.
 */

#include <cmath>
#include <cstring>
#include <iostream>

//#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "util.hpp"

#include "blas.h"

#include "Psp_U.hpp"
#include "compressed_io.hpp"

namespace pwdft {



/* Constructors */

/*******************************************
 *                                         *
 *             Psp_U::Psp_U                *
 *                                         *
 *******************************************/
Psp_U::Psp_U(Ion *myionin, Pneb *mypnebin,
             Strfac *mystrfacin, Control2 &control,
             std::ostream &coutput) 
{
   myion = myionin;
   mypneb = mypnebin;
   mystrfac = mystrfacin;

   uterm   = control.psputerm();
   nuterms = control.pspnuterms();

   if (uterm)
   {

      int nion = myion->nion;
      uterm_l      = new int[nuterms];
      uterm_uscale = new double[nuterms];
      uterm_jscale = new double[nuterms];
      uterm_ions   = new bool[nuterms*nion];
      uterm_vstart = new int[nuterms];

      std::string uterm_output = control.set_psputerm(nion,nuterms,
                                                     uterm_l,
                                                     uterm_uscale,
                                                     uterm_jscale,
                                                     uterm_ions);

      if (mypneb->d3db::parall->is_master())
         coutput << uterm_output << std::endl;

     // for (auto i=0; i<nuterms; ++i)
     // {
     //        call psputerm_gen_vmmmm(int_mb(psputerm_l(1)+i-1),
     //                               2*int_mb(psputerm_l(1)+i-1)+1,
     //                               dbl_mb(psputerm_uscale(1)+i-1),
     //                               dbl_mb(psputerm_jscale(1)+i-1),
     //                               dbl_mb(psputerm_vmmmm(1)
     //                               +int_mb(psputerm_vstart(1)+i-1)))
     // }

   }

}
  

/*******************************************
 *                                         *
 *             Psp_U::v_nonlocal           *
 *                                         *
 *******************************************/

/*******************************************
 *                                         *
 *             Psp_U::f_nonlocal           *
 *                                         *
 *******************************************/

/*******************************************
 *                                         *
 *             Psp_U::gen_vmmmm            *
 *                                         *
 *******************************************/

/*******************************************
 *                                         *
 *             Psp_U::gen_vmm              *
 *                                         *
 *******************************************/

/*******************************************
 *                                         *
 *             Psp_U::gen_l_density        *
 *                                         *
 *******************************************/

/*******************************************
 *                                         *
 *             Psp_U::add_upotential       *
 *                                         *
 *******************************************/

/*******************************************
 *                                         *
 *             Psp_U::uterm_energy         *
 *                                         *
 *******************************************/


} // namespace pwdft

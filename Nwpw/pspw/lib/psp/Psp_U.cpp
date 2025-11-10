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
#include "Paw_gaunt.hpp"
#include "compressed_io.hpp"

namespace pwdft {


/****************************************************
 *                                                  *
 *                uterm_gen_vmmmm                *
 *                                                  *
 ****************************************************/
static void uterm_gen_vmmmm(int l,int lm, double U, double J, double *vmmmm)
{
   double fourpi = 16.0*atan(1.0);
   double F[4];

   if (l==0) 
   {
      F[0] = U;
      F[1] = 0.0;
      F[2] = 0.0;
      F[3] = 0.0;
   }
   if (l==1)
   {
      F[0] = U;
      F[1] = 5*J;
      F[2] = 0.0;
      F[3] = 0.0;
   }
   if (l==2)
   {
      F[0] = U;
      F[1] = 112.0/13.0*J;
      F[2] =  70.0/13.0*J;
      F[3] = 0.0;
   }
   if (l==3)
   {
      F[0] = U;
      F[1] = 4125.0/346.0*J;
      F[2] = 5511.0/692.0*J;
      F[3] = 8151.0/13840*J;
   }

  // Helper to flatten (i,j,k,l) into a linear index
    auto idx4 = [lm](int i, int j, int k, int l4) -> std::size_t {
        // i,j,k,l4 are 0-based within [0..lm-1]
        return static_cast<std::size_t>(((l4 * lm + k) * lm + j) * lm + i);
    };

   // call nwpw_gaunt_init(.false.,3)
   for (auto m1=-l; m1<=l; ++m1)
   for (auto m2=-l; m2<=l; ++m2)
   for (auto m3=-l; m3<=l; ++m3)
   for (auto m4=-l; m4<=l; ++m4)
   {
      const int i1 = m1 + l;
      const int i2 = m2 + l;
      const int i3 = m3 + l;
      const int i4 = m4 + l;

      double value = 0.0;

      Paw_gaunt gnt(true, 3);
      // kk = 1..4  -> k = 0,2,4,6
      for (int kk = 1; kk <= 4; ++kk) 
      {
         const int k = 2 * (kk - 1);
         double aterm = 0.0;

         for (int q = -k; q <= k; ++q) 
         {
            const double g1 = gnt.gaunt(k, q, l, m1, l, m2);
            const double g2 = gnt.gaunt(k, q, l, m4, l, m3);
            aterm += g1 * g2;
         }
         aterm *= fourpi / static_cast<double>(2 * k + 1);

         value += F[kk - 1] * aterm;
      }
 
      vmmmm[idx4(i1,i2,i3,i4)] = value;
   } 

      
}

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

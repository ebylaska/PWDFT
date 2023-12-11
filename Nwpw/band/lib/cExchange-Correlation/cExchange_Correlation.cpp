/*
 *
 *   gga = -1 off
 *
 *   **** LDAs ****
 *   gga = 0 vosko,
 *   gga = 1-9 (reserved for other lda's)
 *
 *   **** GGAs ****
 *   gga = 10 pbe96, pbe
 *   gga = 11 blyp
 *   gga = 12 revpbe
 *   gga = 13 pbesol
 *   gga = 14 hser remainder   (-0.25*Ex(w,pbe,sr) + Ex(pbe) + Ec(pbe))
 *   gga = 15 b3lypr remainder
 *   gga = 16 BEEF
 *   gga = 17 XBEEF-CPBE
 *   gga = 17-99 (reserved for other gga's)
 *
 *   **** hybrids ****
 *   gga = 100-109 (reserved for lda hybrids)
 *   gga = 110  pbe0
 *   gga = 111  ?????
 *   gga = 112  revpbe0
 *   gga = 113  bnl
 *   gga = 114  hse
 *   gga = 115  b3lyp
 *   gga = 116-199 (reserved for hybrids)
 *   gga = 200 hartree-fock
 *
 *   **** meta ggas and hybrid metagga ****
 *   gga = 300 vs98
 *   gga = 301 tpss03
 *   gga = 302 scan
 *   gga = 303 pkzb
 *   gga = 304 m06-l
 *   gga = 305  m06
 *   gga = 306  m06-2x
 */

#include "cExchange_Correlation.hpp"
#include "v_cwexc.hpp"
#include "v_exc.hpp"
#include <algorithm>
#include "parsestring.hpp"

namespace pwdft {

/* Constructors */

/*******************************************
 *                                         *
 *        cXC_Operator::cXC_Operator       *
 *                                         *
 *******************************************/

cXC_Operator::cXC_Operator(Cneb *mygrid, Control2 &control) 
{
   mycneb = mygrid;
   xc_name = control.xc_name();
 
   gga = 0;
 
   if (mystring_contains(mystring_lowercase(xc_name), "vosko")) gga = 0;
   if (mystring_contains(mystring_lowercase(xc_name), "lda"))   gga = 0;

   if (mystring_contains(mystring_lowercase(xc_name), "pbe"))        gga = 10;
   if (mystring_contains(mystring_lowercase(xc_name), "pbe96"))      gga = 10;
   if (mystring_contains(mystring_lowercase(xc_name), "blyp"))       gga = 11;
   if (mystring_contains(mystring_lowercase(xc_name), "revpbe"))     gga = 12;
   if (mystring_contains(mystring_lowercase(xc_name), "pbesol"))     gga = 13;
   if (mystring_contains(mystring_lowercase(xc_name), "hser"))       gga = 14;
   if (mystring_contains(mystring_lowercase(xc_name), "b3lypr"))     gga = 15;
   if (mystring_contains(mystring_lowercase(xc_name), "beef"))       gga = 16;
   if (mystring_contains(mystring_lowercase(xc_name), "xbeef-cpbe")) gga = 17;
   
   if (mystring_contains(mystring_lowercase(xc_name), "pbe0"))    gga = 110;
   if (mystring_contains(mystring_lowercase(xc_name), "blyp0"))   gga = 111;
   if (mystring_contains(mystring_lowercase(xc_name), "revpbe0")) gga = 112;
   if (mystring_contains(mystring_lowercase(xc_name), "bnl"))     gga = 113;
   if (mystring_contains(mystring_lowercase(xc_name), "hse"))     gga = 114;
   if (mystring_contains(mystring_lowercase(xc_name), "b3lyp"))   gga = 115;

   if (mystring_contains(mystring_lowercase(xc_name), "hartree-fock")) gga = 200;
   if (mystring_contains(mystring_lowercase(xc_name), "hf"))           gga = 200;

   use_lda = false;
   use_gga = false;
   use_mgga = false;
   if (gga == 0) {
     use_lda = true;
     xtmp = new double[mycneb->ispin * mycneb->nfft3d];
   }
   if ((gga >= 10) && (gga < 100)) {
     use_gga = true;
     if (mycneb->ispin == 1) {
       rho = new double[mycneb->nfft3d];
 
       grx = new double[mycneb->n2ft3d];
       gry = new double[mycneb->n2ft3d];
       grz = new double[mycneb->n2ft3d];
 
       agr = new double[mycneb->n2ft3d];
       fn = new double[mycneb->n2ft3d];
       fdn = new double[mycneb->n2ft3d];
     } else {
       rho = new double[2 * mycneb->nfft3d];
 
       grx = new double[3 * mycneb->n2ft3d];
       gry = new double[3 * mycneb->n2ft3d];
       grz = new double[3 * mycneb->n2ft3d];
 
       agr = new double[3 * mycneb->n2ft3d];
       fn = new double[2 * mycneb->n2ft3d];
       fdn = new double[3 * mycneb->n2ft3d];
     }
   }
   if ((gga >= 300))
     use_mgga = true;
 
   // std::cout << "xc_name =" << xc_name << std::endl;
}

/*******************************************
 *                                         *
 *        cXC_Operator::v_exc_all           *
 *                                         *
 *******************************************/
void cXC_Operator::v_exc_all(int ispin, double *dn, double *xcp, double *xce) {
  if (use_lda) {
    v_exc(ispin, mycneb->nfft3d, dn, xcp, xce, xtmp);
  } else if (use_gga) {
    v_cwexc(gga, mycneb, dn, 1.0, 1.0, xcp, xce, rho, grx, gry, grz, agr, fn,
            fdn);
  } else if (use_mgga) {
  }
}

} // namespace pwdft

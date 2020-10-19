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


#include        "exchange_correlation.hpp"

/* Constructors */

/*******************************************
 *                                         *
 *        XC_Operator::XC_Operator         *
 *                                         *
 *******************************************/

XC_Operator::XC_Operator(Pneb *mygrid, Control2& control)
{
   mypneb = mygrid;
   xc_name = control.xc_name();

   std::cout << "xc_name =" << xc_name << std::endl;
}



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
 *
 *   **** control dispersion block ****
 *   has_disp = false
 *   is_grimme2 = false
 *   has_vdw  = false
 *   is_vdw2 = false
 *   options_disp = ''
 *  
 *   std::string xc_name,options_disp;
 *   int gga;
 *   bool use_lda, use_gga, use_mgga;
 *
 *   bool has_disp = false;
 *   bool has_vdw  = false;
 *   bool is_grimme2,is_vdw2;
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
 
   // set the grimme and vdw options
   if      (mystring_contains(mystring_lowercase(xc_name), "-grimme2")) { has_disp = true; is_grimme2 = true;  options_disp = " -old -noprint";}
   else if (mystring_contains(mystring_lowercase(xc_name), "-grimme3")) { has_disp = true; is_grimme2 = false; options_disp = " -zero -noprint";}
   else if (mystring_contains(mystring_lowercase(xc_name), "-grimme4")) { has_disp = true; is_grimme2 = false; options_disp = " -bj -num -noprint";}
   else if (mystring_contains(mystring_lowercase(xc_name), "-grimme5")) { has_disp = true; is_grimme2 = false; options_disp = " -zerom -noprint";}
   else if (mystring_contains(mystring_lowercase(xc_name), "-grimme6")) { has_disp = true; is_grimme2 = false; options_disp = " -bjm -num -noprint";}
   else if (mystring_contains(mystring_lowercase(xc_name), "-vdw2"))    { has_vdw = true; is_vdw2 = false; }
   else if (mystring_contains(mystring_lowercase(xc_name), "-vdw"))     { has_vdw = true; is_vdw2 = true; }

   gga = 0;
 
   if (mystring_contains(mystring_lowercase(xc_name), "vosko")) gga = 0;
   if (mystring_contains(mystring_lowercase(xc_name), "lda"))   gga = 0;

   // set the gga options
   if (mystring_contains(mystring_lowercase(xc_name), "pbe"))        {gga = 10; if (has_disp) options_disp = "-func pbe" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "pbe96"))      {gga = 10; if (has_disp) options_disp = "-func pbe" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "blyp"))       {gga = 11; if (has_disp) options_disp = "-func b-lyp" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "revpbe"))     {gga = 12; if (has_disp) options_disp = "-func revpbe" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "pbesol"))     {gga = 13; if (has_disp) options_disp = "-func pbesol" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "hser"))       {gga = 14; if (has_disp) options_disp = "-func hse06" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "b3lypr"))     {gga = 15; if (has_disp) options_disp = "-func b3-lyp" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "beef"))       {gga = 16; if (has_disp) options_disp = "-func pbesol" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "xbeef-cpbe")) {gga = 17; if (has_disp) options_disp = "-func pbesol" + options_disp;}

   if (mystring_contains(mystring_lowercase(xc_name), "pbe0"))    {gga = 110; if (has_disp) options_disp = "-func pbe0" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "blyp0"))   {gga = 111; if (has_disp) options_disp = "-func b3-lyp" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "revpbe0")) {gga = 112; if (has_disp) options_disp = "-func revpbe0" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "bnl"))     {gga = 113; if (has_disp) options_disp = "-func hse06" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "hse"))     {gga = 114; if (has_disp) options_disp = "-func hse06" + options_disp;}
   if (mystring_contains(mystring_lowercase(xc_name), "b3lyp"))   {gga = 115; if (has_disp) options_disp = "-func b3-lyp" + options_disp;}

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
       rho = new double[mycneb->n2ft3d]; // real
 
       grx = new double[mycneb->n2ft3d]; // complex
       gry = new double[mycneb->n2ft3d]; // complex
       grz = new double[mycneb->n2ft3d]; // complex
 
       agr = new double[mycneb->n2ft3d]; // real|complex
       fn  = new double[mycneb->n2ft3d]; // real|complex
       fdn = new double[mycneb->n2ft3d]; // real|complex
     } else {
       rho = new double[2 * mycneb->n2ft3d]; // real
 
       grx = new double[3 * mycneb->n2ft3d]; // complex
       gry = new double[3 * mycneb->n2ft3d]; // complex
       grz = new double[3 * mycneb->n2ft3d]; // complex
 
       agr = new double[3 * mycneb->n2ft3d]; // real|complex
       fn  = new double[2 * mycneb->n2ft3d]; // real|complex
       fdn = new double[3 * mycneb->n2ft3d]; // real|complex
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
    //std::cout << "dn=" << dn[0] << " " << dn[1] << std::endl;
    //std::cout << "xcp=" << xcp[0] << " " << xcp[1] << std::endl;
    //double sumall = 0.0;
    //for (auto i=0; i<mycneb->nfft3d; ++i)
   // {
   //    std::cout << "i=" << i << " dnall=" << dn[i] << " xcp=" << xcp[i] << std::endl;
   //    sumall += dn[i];
   // }
   // std::cout << "sumall=" << sumall << std::endl;

  } else if (use_gga) {
    v_cwexc(gga, mycneb, dn, 1.0, 1.0, xcp, xce, rho, grx, gry, grz, agr, fn,
            fdn);
  } else if (use_mgga) {
  }
}

} // namespace pwdft

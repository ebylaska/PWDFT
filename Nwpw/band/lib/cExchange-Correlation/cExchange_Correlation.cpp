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

   const std::string name = mystring_lowercase(xc_name);

   auto add_disp = [&](const std::string& func_flag) {
      if (!has_disp) return;
      // ensure exactly one -func … prefix
      if (!func_flag.empty()) {
         // always keep a leading space before subsequent flags
         options_disp = "-func " + func_flag + (options_disp.empty() ? "" : " " + options_disp);
      }
   };


   gga = 0;

   // set the grimme and vdw options
   if      (mystring_contains(name, "-grimme2")) {has_disp = true; is_grimme2 = true;  options_disp = "-old -noprint";}
   else if (mystring_contains(name, "-grimme3")) {has_disp = true; is_grimme2 = false; options_disp = "-zero -noprint";}
   else if (mystring_contains(name, "-grimme4")) {has_disp = true; is_grimme2 = false; options_disp = "-bj -num -noprint";}
   else if (mystring_contains(name, "-grimme5")) {has_disp = true; is_grimme2 = false; options_disp = "-zerom -noprint";}
   else if (mystring_contains(name, "-grimme6")) {has_disp = true; is_grimme2 = false; options_disp = "-bjm -num -noprint";}

   if      (mystring_contains(name, "-vdw2"))    {has_vdw = true; is_vdw2 = true; }
   else if (mystring_contains(name, "-vdw"))     {has_vdw = true; is_vdw2 = false; }

   // set the gga options

   // ---- XC family: choose longest/specific first; first match wins ----
   // Hybrids first (specific → general)
   if      (mystring_contains(name, "revpbe0")) {gga = 112; add_disp("revpbe0");}
   else if (mystring_contains(name, "pbe0"))    {gga = 110; add_disp("pbe0");}
   else if (mystring_contains(name, "hse"))     {gga = 114; add_disp("hse06");}
   else if (mystring_contains(name, "bnl"))     {gga = 113; add_disp("hse06");}  // if that’s really what you want
   else if (mystring_contains(name, "b3lypr"))  {gga = 115; add_disp("b3-lyp");} // treat b3lypr as hybrid remainder path
   else if (mystring_contains(name, "blyp0"))   {gga = 111; add_disp("b3-lyp");} // if this alias is desired

   // Non-hybrid GGAs (specific → general)
   else if (mystring_contains(name, "revpbe"))     {gga = 12; add_disp("revpbe");}
   else if (mystring_contains(name, "pbesol"))     {gga = 13; add_disp("pbesol");}
   else if (mystring_contains(name, "pbe96"))      {gga = 10; add_disp("pbe");}
   else if (mystring_contains(name, "xbeef-cpbe")) {gga = 17; add_disp("pbesol");} // confirm
   else if (mystring_contains(name, "beef"))       {gga = 16; add_disp("pbesol");} // confirm
   else if (mystring_contains(name, "blyp"))       {gga = 11; add_disp("b-lyp");}
   else if (mystring_contains(name, "pbe"))        {gga = 10; add_disp("pbe");}

   // LDA family
   else if (mystring_contains(name, "vosko") || mystring_contains(name, "lda")) {gga = 0;}

   // HF exact exchange - dispersion func flag generally not needed/used here
   if (mystring_contains(name, "hartree-fock") || name == "hf" || name.find(" hf ") != std::string::npos) {gga = 200;}

   // Meta-GGAs - future options
   if      (mystring_contains(name, "scan"))   gga = 302;
   else if (mystring_contains(name, "tpss03")) gga = 301;
   else if (mystring_contains(name, "vs98"))   gga = 300;
   else if (mystring_contains(name, "pkzb"))   gga = 303;
   else if (mystring_contains(name, "m06-2x")) gga = 306;
   else if (mystring_contains(name, "m06-l"))  gga = 304;
   else if (mystring_contains(name, "m06"))    gga = 305;


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

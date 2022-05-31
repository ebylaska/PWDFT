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


#include	<algorithm>
#include	"v_exc.hpp"
#include	"v_bwexc.hpp"
#include	"exchange_correlation.hpp"

namespace pwdft {


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

   std::transform(xc_name.begin(), xc_name.end(), xc_name.begin(), ::tolower);

   gga = 0;

   if (xc_name.compare("vosko") == 0) gga=0;
   if (xc_name.compare("lda")   == 0) gga=0;

   if (xc_name.compare("pbe")        == 0) gga=10;
   if (xc_name.compare("pbe96")      == 0) gga=10;
   if (xc_name.compare("blyp")       == 0) gga=11;
   if (xc_name.compare("revpbe")     == 0) gga=12;
   if (xc_name.compare("pbesol")     == 0) gga=13;
   if (xc_name.compare("hser")       == 0) gga=14;
   if (xc_name.compare("b3lypr")     == 0) gga=15;
   if (xc_name.compare("beef")       == 0) gga=16;
   if (xc_name.compare("xbeef-cpbe") == 0) gga=17;

   if (xc_name.compare("pbe0")    == 0) gga=110;
   if (xc_name.compare("blyp0")   == 0) gga=111;
   if (xc_name.compare("revpbe0") == 0) gga=112;
   if (xc_name.compare("bnl")     == 0) gga=113;
   if (xc_name.compare("hse")     == 0) gga=114;
   if (xc_name.compare("b3lyp")   == 0) gga=115;

   if (xc_name.compare("hartree-fock") == 0) gga=200;
   if (xc_name.compare("hf")           == 0) gga=200;

   use_lda  = false;
   use_gga  = false;
   use_mgga = false;
   if (gga==0) {
      use_lda = true;
      xtmp = new double [mypneb->ispin*mypneb->n2ft3d];
   }
   if ((gga>=10) && (gga<100)) {
        use_gga  = true;
        if (mypneb->ispin==1) {
           rho = new double [mypneb->n2ft3d];

           grx = new double [mypneb->n2ft3d];
           gry = new double [mypneb->n2ft3d];
           grz = new double [mypneb->n2ft3d];

           agr = new double [mypneb->n2ft3d];
           fn  = new double [mypneb->n2ft3d];
           fdn = new double [mypneb->n2ft3d];
        } else {
           rho = new double [2*mypneb->n2ft3d];

           grx = new double [3*mypneb->n2ft3d];
           gry = new double [3*mypneb->n2ft3d];
           grz = new double [3*mypneb->n2ft3d];

           agr = new double [3*mypneb->n2ft3d];
           fn  = new double [2*mypneb->n2ft3d];
           fdn = new double [3*mypneb->n2ft3d];
        }
   }
   if ((gga>=300))             use_mgga = true;

    
   //std::cout << "xc_name =" << xc_name << std::endl;
}


/*******************************************
 *                                         *
 *        XC_Operator::v_exc_all           *
 *                                         *
 *******************************************/
void XC_Operator::v_exc_all(int ispin, double *dn, double *xcp, double *xce) {
   if (use_lda) {
      v_exc(ispin,mypneb->n2ft3d,dn,xcp,xce,xtmp);
   } else if (use_gga) {
      v_bwexc(gga,mypneb,dn,1.0,1.0,
              xcp,xce,
              rho,grx,gry,grz,agr,fn,fdn);
   } else if (use_mgga) {
   }
}


}

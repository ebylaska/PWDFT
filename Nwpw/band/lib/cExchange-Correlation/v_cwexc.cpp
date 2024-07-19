
#include "Cneb.hpp"
#include "b3lyp.hpp"
#include "beef_gga.hpp"
#include "blyp.hpp"
#include "hsepbe.hpp"
#include "pbe96.hpp"
#include "pbesol.hpp"
#include "revpbe.hpp"

namespace pwdft {

#define dncut 1.0e-30

/********************************
 *				                    *
 *            v_cwexc           *
 *                              *
 ********************************/
// rho,agr,fn,fdn - real
// grx,gry,grz - complex
void v_cwexc(const int gga, Cneb *mycneb, const double *dn,
             const double x_parameter, const double c_parameter, double *xcp,
             double *xce, double *rho, double *grx, double *gry, double *grz,
             double *agr, double *fn, double *fdn) 
{
   double *rhog = fn; // complex
   double *Gx = mycneb->Gpackxyz(0, 0);
   double *Gy = mycneb->Gpackxyz(0, 1);
   double *Gz = mycneb->Gpackxyz(0, 2);
   double scal1 = 1.0 / ((double)((mycneb->nx) * (mycneb->ny) * (mycneb->nz)));
 
   /**********************************
    ***** restricted calculation *****
    **********************************/
   if (mycneb->ispin == 1) 
   {
      // mycneb->r_nzero(2,agr);
      mycneb->r_zero(agr);
     
      /* calculate rho tmp1=rho(g) */
      mycneb->rr_copy(dn,rho);
      mycneb->rr_Sum(dn,rho);
      mycneb->rc_SMul(scal1,rho,rhog);
      mycneb->rc_fft3d(rhog);
      mycneb->c_pack(0,rhog);
     
      /* calculate gr = grad n - all complex*/
      mycneb->tcc_pack_iMul(0, Gx, rhog, grx);
      mycneb->tcc_pack_iMul(0, Gy, rhog, gry);
      mycneb->tcc_pack_iMul(0, Gz, rhog, grz);
      mycneb->c_unpack(0, grx);
      mycneb->c_unpack(0, gry);
      mycneb->c_unpack(0, grz);
      mycneb->cr_fft3d(grx);
      mycneb->cr_fft3d(gry);
      mycneb->cr_fft3d(grz);
     
      /* calculate agr = |grad n| -- agr needs to be real */
      mycneb->cr_sqr(grx, agr);
      mycneb->cr_addsqr(gry, agr);
      mycneb->cr_addsqr(grz, agr);
      mycneb->r_sqrt(agr);
     
      switch (gga) {
      case 10:
        gen_PBE96_BW_restricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 11:
        gen_BLYP_BW_restricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 12:
        gen_revPBE_BW_restricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 13:
        gen_PBEsol_BW_restricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 14:
        gen_HSE_BW_restricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 15:
        gen_B3LYP_BW_restricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 16:
        gen_BEEF_BW_restricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, 0.6001664769, xce, fn, fdn);
        break;
      case 17:
        gen_BEEF_BW_restricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, 0.0, xce, fn, fdn);
        break;
     
      default:
        gen_PBE96_BW_restricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
      }
     
      /* calculate df/d|grad n| *(grad n)/|grad n| */
      mycneb->rc_Divide(agr, grx);
      mycneb->rc_Divide(agr, gry);
      mycneb->rc_Divide(agr, grz);

      mycneb->rc_Mul(fdn, grx);
      mycneb->rc_Mul(fdn, gry);
      mycneb->rc_Mul(fdn, grz);
     
      mycneb->c_SMul(scal1, grx);
      mycneb->c_SMul(scal1, gry);
      mycneb->c_SMul(scal1, grz);
     
      mycneb->rc_fft3d(grx);
      mycneb->rc_fft3d(gry);
      mycneb->rc_fft3d(grz);
     
      mycneb->c_pack(0, grx);
      mycneb->c_pack(0, gry);
      mycneb->c_pack(0, grz);
     
      mycneb->rc_pack_iMul(0, Gx, grx);
      mycneb->rc_pack_iMul(0, Gy, gry);
      mycneb->rc_pack_iMul(0, Gz, grz);
     
      mycneb->cccc_pack_Sum(0, grx, gry, grz, fdn);
      mycneb->c_unpack(0, fdn);
      mycneb->cr_fft3d(fdn);

      mycneb->rcr_Minus(fn, fdn, xcp);
   }
 
   /************************************
    ***** unrestricted calculation *****
    ************************************/
   else 
   {
      mycneb->r_nzero(3, agr);
     
      double *rhoup = rho;
      double *rhodn = rho + mycneb->nfft3d;
      double *grupx = grx;
      double *grupy = gry;
      double *grupz = grz;
      double *grdnx = grx + mycneb->n2ft3d;
      double *grdny = gry + mycneb->n2ft3d;
      double *grdnz = grz + mycneb->n2ft3d;
      double *grallx = grx + 2 * mycneb->n2ft3d;
      double *grally = gry + 2 * mycneb->n2ft3d;
      double *grallz = grz + 2 * mycneb->n2ft3d;
      double *agrup = agr;
      double *agrdn = agr + mycneb->nfft3d;
      double *agrall = agr + 2 * mycneb->nfft3d;
     
      // double *dnup = dn;
      // double *dndn = dn + mycneb->n2ft3d;
     
      double *xcpup = xcp;
      double *xcpdn = xcp + mycneb->nfft3d;
     
      double *fnup = fn;
      double *fndn = fn + mycneb->nfft3d;
      ;
     
      double *fdnup = fdn;
      double *fdndn = fdn + mycneb->nfft3d;
      ;
      double *fdnall = fdn + 2 * mycneb->nfft3d;
      ;
     
      /* calculate rhoup  */
      mycneb->rr_copy(dn, rhoup);
      mycneb->rc_SMul(scal1, rhoup, rhog);
      mycneb->rc_fft3d(rhog);
      mycneb->c_pack(0, rhog);
     
      /* calculate   grup= grad nup */
      mycneb->tcc_pack_iMul(0, Gx, rhog, grupx);
      mycneb->tcc_pack_iMul(0, Gy, rhog, grupy);
      mycneb->tcc_pack_iMul(0, Gz, rhog, grupz);
      mycneb->c_unpack(0, grupx);
      mycneb->c_unpack(0, grupy);
      mycneb->c_unpack(0, grupz);
      mycneb->cr_fft3d(grupx);
      mycneb->cr_fft3d(grupy);
      mycneb->cr_fft3d(grupz);
     
      /* calculate agrup = |grad nup| */
      mycneb->cr_sqr(grupx, agrup);
      mycneb->cr_addsqr(grupy, agrup);
      mycneb->cr_addsqr(grupz, agrup);
      mycneb->r_sqrt(agrup);
     
      /* calculate rhodn  */
      mycneb->rr_copy(dn + mycneb->nfft3d, rhodn);
      mycneb->rc_SMul(scal1, rhodn, rhog);
      mycneb->rc_fft3d(rhog);
      mycneb->c_pack(0, rhog);
     
      /* calculate   grdn= grad ndn */
      mycneb->tcc_pack_iMul(0, Gx, rhog, grdnx);
      mycneb->tcc_pack_iMul(0, Gy, rhog, grdny);
      mycneb->tcc_pack_iMul(0, Gz, rhog, grdnz);
      mycneb->c_unpack(0, grdnx);
      mycneb->c_unpack(0, grdny);
      mycneb->c_unpack(0, grdnz);
      mycneb->cr_fft3d(grdnx);
      mycneb->cr_fft3d(grdny);
      mycneb->cr_fft3d(grdnz);
     
      /* calculate agrdn = |grad ndn| */
      mycneb->cr_sqr(grdnx, agrdn);
      mycneb->cr_addsqr(grdny, agrdn);
      mycneb->cr_addsqr(grdnz, agrdn);
      mycneb->r_sqrt(agrdn);
     
      /* calculate agrall = |grad nup +grad ndn| */
      mycneb->ccc_Sum(grupx, grdnx, grallx);
      mycneb->ccc_Sum(grupy, grdny, grally);
      mycneb->ccc_Sum(grupz, grdnz, grallz);
      mycneb->cr_sqr(grallx, agrall);
      mycneb->cr_addsqr(grally, agrall);
      mycneb->cr_addsqr(grallz, agrall);
      mycneb->r_sqrt(agrall);
     
      switch (gga) {
      case 10:
        gen_PBE96_BW_unrestricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 11:
        gen_BLYP_BW_unrestricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 12:
        gen_revPBE_BW_unrestricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 13:
        gen_PBEsol_BW_unrestricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 14:
        gen_HSE_BW_unrestricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 15:
        gen_B3LYP_BW_unrestricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 16:
        gen_BEEF_BW_unrestricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, 0.6001664769, xce, fn, fdn);
        break;
      case 17:
        gen_BEEF_BW_unrestricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, 0.0, xce, fn, fdn);
        break;
     
      default:
        gen_PBE96_BW_unrestricted(mycneb->nfft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
      }
     
      /**** calculate df/d|grad nup|* (grad nup)/|grad nup|  ****
       **** calculate df/d|grad ndn|* (grad ndn)/|grad ndn|  ****
       **** calculate df/d|grad n|  * (grad n)/|grad n|  ****/
      mycneb->rc_Divide(agrup, grupx);
      mycneb->rc_Divide(agrup, grupy);
      mycneb->rc_Divide(agrup, grupz);
      mycneb->rc_Divide(agrdn, grdnx);
      mycneb->rc_Divide(agrdn, grdny);
      mycneb->rc_Divide(agrdn, grdnz);
      mycneb->rc_Divide(agrall, grallx);
      mycneb->rc_Divide(agrall, grally);
      mycneb->rc_Divide(agrall, grallz);
     
      mycneb->rc_Mul(fdnup, grupx);
      mycneb->rc_Mul(fdnup, grupy);
      mycneb->rc_Mul(fdnup, grupz);
      mycneb->rc_Mul(fdndn, grdnx);
      mycneb->rc_Mul(fdndn, grdny);
      mycneb->rc_Mul(fdndn, grdnz);
      mycneb->rc_Mul(fdnall, grallx);
      mycneb->rc_Mul(fdnall, grally);
      mycneb->rc_Mul(fdnall, grallz);
     
      /**** calculate (df/d|grad nup|* (grad nup)/|grad nup|)  ****
       ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****
       **** calculate (df/d|grad ndn|* (grad ndn)/|grad ndn|)  ****
       ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****/
      mycneb->cc_Sum(grallx, grupx);
      mycneb->cc_Sum(grally, grupy);
      mycneb->cc_Sum(grallz, grupz);
      mycneb->cc_Sum(grallx, grdnx);
      mycneb->cc_Sum(grally, grdny);
      mycneb->cc_Sum(grallz, grdnz);
     
      mycneb->c_SMul(scal1, grupx);
      mycneb->c_SMul(scal1, grupy);
      mycneb->c_SMul(scal1, grupz);
      mycneb->c_SMul(scal1, grdnx);
      mycneb->c_SMul(scal1, grdny);
      mycneb->c_SMul(scal1, grdnz);
     
      /* put sums by G-space */
      mycneb->rc_fft3d(grupx);
      mycneb->rc_fft3d(grupy);
      mycneb->rc_fft3d(grupz);
      mycneb->rc_fft3d(grdnx);
      mycneb->rc_fft3d(grdny);
      mycneb->rc_fft3d(grdnz);
     
      mycneb->c_pack(0, grupx);
      mycneb->c_pack(0, grupy);
      mycneb->c_pack(0, grupz);
      mycneb->c_pack(0, grdnx);
      mycneb->c_pack(0, grdny);
      mycneb->c_pack(0, grdnz);
     
      /* multiply sums by G vector */
      mycneb->rc_pack_iMul(0, Gx, grupx);
      mycneb->rc_pack_iMul(0, Gy, grupy);
      mycneb->rc_pack_iMul(0, Gz, grupz);
      mycneb->rc_pack_iMul(0, Gx, grdnx);
      mycneb->rc_pack_iMul(0, Gy, grdny);
      mycneb->rc_pack_iMul(0, Gz, grdnz);
     
      /* addup dot products */
      mycneb->cccc_pack_Sum(0, grupx, grupy, grupz, fdnup);
      mycneb->cccc_pack_Sum(0, grdnx, grdny, grdnz, fdndn);
     
      /* put back in r-space and subtract from df/dnup,df/dndn */
      mycneb->c_unpack(0, fdnup);
      mycneb->c_unpack(0, fdndn);
      mycneb->cr_fft3d(fdnup);
      mycneb->cr_fft3d(fdndn);
      mycneb->rcr_Minus(fnup, fdnup, xcpup);
      mycneb->rcr_Minus(fndn, fdndn, xcpdn);
   }
}

} // namespace pwdft

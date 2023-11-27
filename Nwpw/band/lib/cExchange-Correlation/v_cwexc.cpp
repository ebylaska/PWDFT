
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
void v_cwexc(const int gga, Cneb *mycneb, const double *dn,
             const double x_parameter, const double c_parameter, double *xcp,
             double *xce, double *rho, double *grx, double *gry, double *grz,
             double *agr, double *fn, double *fdn) 
{
   double *rhog = fn;
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
      mycneb->rr_copy(dn, rho);
      mycneb->rr_Sum(dn, rho);
      mycneb->rr_SMul(scal1, rho, rhog);
      mycneb->rc_fft3d(rhog);
      mycneb->c_pack(0, rhog);
     
      /* calculate gr = grad n */
      mycneb->tcc_pack_iMul(0, Gx, rhog, grx);
      mycneb->tcc_pack_iMul(0, Gy, rhog, gry);
      mycneb->tcc_pack_iMul(0, Gz, rhog, grz);
      mycneb->c_unpack(0, grx);
      mycneb->c_unpack(0, gry);
      mycneb->c_unpack(0, grz);
      mycneb->cr_fft3d(grx);
      mycneb->cr_fft3d(gry);
      mycneb->cr_fft3d(grz);
     
      /* calculate agr = |grad n| */
      mycneb->rr_sqr(grx, agr);
      mycneb->rr_addsqr(gry, agr);
      mycneb->rr_addsqr(grz, agr);
      mycneb->r_sqrt(agr);
     
      switch (gga) {
      case 10:
        gen_PBE96_BW_restricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
        break;
      case 11:
        gen_BLYP_BW_restricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                               xce, fn, fdn);
        break;
      case 12:
        gen_revPBE_BW_restricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, xce, fn, fdn);
        break;
      case 13:
        gen_PBEsol_BW_restricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, xce, fn, fdn);
        break;
      case 14:
        gen_HSE_BW_restricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                              xce, fn, fdn);
        break;
      case 15:
        gen_B3LYP_BW_restricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
        break;
      case 16:
        gen_BEEF_BW_restricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                               0.6001664769, xce, fn, fdn);
        break;
      case 17:
        gen_BEEF_BW_restricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                               0.0, xce, fn, fdn);
        break;
     
      default:
        gen_PBE96_BW_restricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
      }
     
      /* calculate df/d|grad n| *(grad n)/|grad n| */
      mycneb->rr_Divide(agr, grx);
      mycneb->rr_Divide(agr, gry);
      mycneb->rr_Divide(agr, grz);
      mycneb->rr_Mul(fdn, grx);
      mycneb->rr_Mul(fdn, gry);
      mycneb->rr_Mul(fdn, grz);
     
      mycneb->r_SMul(scal1, grx);
      mycneb->r_SMul(scal1, gry);
      mycneb->r_SMul(scal1, grz);
     
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
      mycneb->rrr_Minus(fn, fdn, xcp);
   }
 
   /************************************
    ***** unrestricted calculation *****
    ************************************/
   else 
   {
      mycneb->r_nzero(3, agr);
     
      double *rhoup = rho;
      double *rhodn = rho + mycneb->n2ft3d;
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
      double *agrdn = agr + mycneb->n2ft3d;
      double *agrall = agr + 2 * mycneb->n2ft3d;
     
      // double *dnup = dn;
      // double *dndn = dn + mycneb->n2ft3d;
     
      double *xcpup = xcp;
      double *xcpdn = xcp + mycneb->n2ft3d;
     
      double *fnup = fn;
      double *fndn = fn + mycneb->n2ft3d;
      ;
     
      double *fdnup = fdn;
      double *fdndn = fdn + mycneb->n2ft3d;
      ;
      double *fdnall = fdn + 2 * mycneb->n2ft3d;
      ;
     
      /* calculate rhoup  */
      mycneb->rr_copy(dn, rhoup);
      mycneb->rr_SMul(scal1, rhoup, rhog);
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
      mycneb->rr_sqr(grupx, agrup);
      mycneb->rr_addsqr(grupy, agrup);
      mycneb->rr_addsqr(grupz, agrup);
      mycneb->r_sqrt(agrup);
     
      /* calculate rhodn  */
      mycneb->rr_copy(dn + mycneb->n2ft3d, rhodn);
      mycneb->rr_SMul(scal1, rhodn, rhog);
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
      mycneb->rr_sqr(grdnx, agrdn);
      mycneb->rr_addsqr(grdny, agrdn);
      mycneb->rr_addsqr(grdnz, agrdn);
      mycneb->r_sqrt(agrdn);
     
      /* calculate agrall = |grad nup +grad ndn| */
      mycneb->rrr_Sum(grupx, grdnx, grallx);
      mycneb->rrr_Sum(grupy, grdny, grally);
      mycneb->rrr_Sum(grupz, grdnz, grallz);
      mycneb->rr_sqr(grallx, agrall);
      mycneb->rr_addsqr(grally, agrall);
      mycneb->rr_addsqr(grallz, agrall);
      mycneb->r_sqrt(agrall);
     
      switch (gga) {
      case 10:
        gen_PBE96_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                  c_parameter, xce, fn, fdn);
        break;
      case 11:
        gen_BLYP_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, xce, fn, fdn);
        break;
      case 12:
        gen_revPBE_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                   c_parameter, xce, fn, fdn);
        break;
      case 13:
        gen_PBEsol_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                   c_parameter, xce, fn, fdn);
        break;
      case 14:
        gen_HSE_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
        break;
      case 15:
        gen_B3LYP_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                  c_parameter, xce, fn, fdn);
        break;
      case 16:
        gen_BEEF_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, 0.6001664769, xce, fn, fdn);
        break;
      case 17:
        gen_BEEF_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, 0.0, xce, fn, fdn);
        break;
     
      default:
        gen_PBE96_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter,
                                  c_parameter, xce, fn, fdn);
      }
     
      /**** calculate df/d|grad nup|* (grad nup)/|grad nup|  ****
       **** calculate df/d|grad ndn|* (grad ndn)/|grad ndn|  ****
       **** calculate df/d|grad n|  * (grad n)/|grad n|  ****/
      mycneb->rr_Divide(agrup, grupx);
      mycneb->rr_Divide(agrup, grupy);
      mycneb->rr_Divide(agrup, grupz);
      mycneb->rr_Divide(agrdn, grdnx);
      mycneb->rr_Divide(agrdn, grdny);
      mycneb->rr_Divide(agrdn, grdnz);
      mycneb->rr_Divide(agrall, grallx);
      mycneb->rr_Divide(agrall, grally);
      mycneb->rr_Divide(agrall, grallz);
     
      mycneb->rr_Mul(fdnup, grupx);
      mycneb->rr_Mul(fdnup, grupy);
      mycneb->rr_Mul(fdnup, grupz);
      mycneb->rr_Mul(fdndn, grdnx);
      mycneb->rr_Mul(fdndn, grdny);
      mycneb->rr_Mul(fdndn, grdnz);
      mycneb->rr_Mul(fdnall, grallx);
      mycneb->rr_Mul(fdnall, grally);
      mycneb->rr_Mul(fdnall, grallz);
     
      /**** calculate (df/d|grad nup|* (grad nup)/|grad nup|)  ****
       ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****
       **** calculate (df/d|grad ndn|* (grad ndn)/|grad ndn|)  ****
       ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****/
      mycneb->rr_Sum(grallx, grupx);
      mycneb->rr_Sum(grally, grupy);
      mycneb->rr_Sum(grallz, grupz);
      mycneb->rr_Sum(grallx, grdnx);
      mycneb->rr_Sum(grally, grdny);
      mycneb->rr_Sum(grallz, grdnz);
     
      mycneb->r_SMul(scal1, grupx);
      mycneb->r_SMul(scal1, grupy);
      mycneb->r_SMul(scal1, grupz);
      mycneb->r_SMul(scal1, grdnx);
      mycneb->r_SMul(scal1, grdny);
      mycneb->r_SMul(scal1, grdnz);
     
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
      mycneb->rrr_Minus(fnup, fdnup, xcpup);
      mycneb->rrr_Minus(fndn, fdndn, xcpdn);
   }
}

} // namespace pwdft

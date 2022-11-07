      
#include        "Pneb.hpp"
#include        "pbe96.hpp"
#include        "blyp.hpp"
#include        "revpbe.hpp"
#include        "pbesol.hpp"
#include        "beef_gga.hpp"
#include        "hsepbe.hpp"
#include        "b3lyp.hpp"

namespace pwdft {


#define dncut	1.0e-30

/********************************
 *				*
 *            v_bwexc		*
 *				*
 ********************************/
void v_bwexc(const int gga, Pneb *mypneb,
             const double *dn, 
             const double x_parameter, const double c_parameter,
             double *xcp, double *xce,
             double *rho, double *grx, double *gry, double *grz,
             double *agr, double *fn, double *fdn)
{

  double *rhog = fn;
  double *Gx = mypneb->Gpackxyz(0,0);
  double *Gy = mypneb->Gpackxyz(0,1);
  double *Gz = mypneb->Gpackxyz(0,2);
  double scal1 = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));


  /********************************** 
   ***** restricted calculation ***** 
   **********************************/
   if (mypneb->ispin==1)
   {
      //mypneb->r_nzero(2,agr);
      mypneb->r_zero(agr);

      /* calculate rho tmp1=rho(g) */
      mypneb->rr_copy(dn,rho);
      mypneb->rr_Sum(dn,rho);
      mypneb->rr_SMul(scal1,rho,rhog);
      mypneb->rc_fft3d(rhog);
      mypneb->c_pack(0,rhog);

      /* calculate gr = grad n */
      mypneb->tcc_pack_iMul(0,Gx,rhog,grx);
      mypneb->tcc_pack_iMul(0,Gy,rhog,gry);
      mypneb->tcc_pack_iMul(0,Gz,rhog,grz);
      mypneb->c_unpack(0,grx);
      mypneb->c_unpack(0,gry);
      mypneb->c_unpack(0,grz);
      mypneb->cr_fft3d(grx);
      mypneb->cr_fft3d(gry);
      mypneb->cr_fft3d(grz);

      /* calculate agr = |grad n| */
      mypneb->rr_sqr(grx,agr);
      mypneb->rr_addsqr(gry,agr);
      mypneb->rr_addsqr(grz,agr);
      mypneb->r_sqrt(agr);

      switch (gga) {
         case 10 :
            gen_PBE96_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 11 :
            gen_BLYP_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 12 :
            gen_revPBE_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 13 :
            gen_PBEsol_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 14 :
            gen_HSE_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 15 :
            gen_B3LYP_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 16 :
            gen_BEEF_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,0.6001664769,xce,fn,fdn);
            break;
         case 17 :
            gen_BEEF_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,0.0,xce,fn,fdn);
            break;

         default:
            gen_PBE96_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
      }

      /* calculate df/d|grad n| *(grad n)/|grad n| */
      mypneb->rr_Divide(agr,grx);
      mypneb->rr_Divide(agr,gry);
      mypneb->rr_Divide(agr,grz);
      mypneb->rr_Mul(fdn,grx);
      mypneb->rr_Mul(fdn,gry);
      mypneb->rr_Mul(fdn,grz);

      mypneb->r_SMul(scal1,grx);
      mypneb->r_SMul(scal1,gry);
      mypneb->r_SMul(scal1,grz);

      mypneb->rc_fft3d(grx);
      mypneb->rc_fft3d(gry);
      mypneb->rc_fft3d(grz);

      mypneb->c_pack(0,grx);
      mypneb->c_pack(0,gry);
      mypneb->c_pack(0,grz);

      mypneb->tc_pack_iMul(0,Gx,grx);
      mypneb->tc_pack_iMul(0,Gy,gry);
      mypneb->tc_pack_iMul(0,Gz,grz);

      mypneb->cccc_pack_Sum(0,grx,gry,grz,fdn);
      mypneb->c_unpack(0,fdn);
      mypneb->cr_fft3d(fdn);
      mypneb->rrr_Minus(fn,fdn,xcp);
   } 

  /************************************ 
   ***** unrestricted calculation ***** 
   ************************************/
   else 
   {
      mypneb->r_nzero(3,agr);

      double *rhoup = rho;
      double *rhodn = rho + mypneb->n2ft3d;
      double *grupx  = grx;
      double *grupy  = gry;
      double *grupz  = grz;
      double *grdnx  = grx + mypneb->n2ft3d;
      double *grdny  = gry + mypneb->n2ft3d;
      double *grdnz  = grz + mypneb->n2ft3d;
      double *grallx = grx + 2*mypneb->n2ft3d;
      double *grally = gry + 2*mypneb->n2ft3d;
      double *grallz = grz + 2*mypneb->n2ft3d;
      double *agrup  = agr;
      double *agrdn  = agr + mypneb->n2ft3d;
      double *agrall = agr + 2*mypneb->n2ft3d;

      //double *dnup = dn;
      //double *dndn = dn + mypneb->n2ft3d;

      double *xcpup = xcp;
      double *xcpdn = xcp + mypneb->n2ft3d;

      double *fnup = fn;
      double *fndn = fn + mypneb->n2ft3d;;

      double *fdnup  = fdn;
      double *fdndn  = fdn + mypneb->n2ft3d;;
      double *fdnall = fdn + 2*mypneb->n2ft3d;;

      /* calculate rhoup  */
      mypneb->rr_copy(dn,rhoup);
      mypneb->rr_SMul(scal1,rhoup,rhog);
      mypneb->rc_fft3d(rhog);
      mypneb->c_pack(0,rhog);

      /* calculate   grup= grad nup */
      mypneb->tcc_pack_iMul(0,Gx,rhog,grupx);
      mypneb->tcc_pack_iMul(0,Gy,rhog,grupy);
      mypneb->tcc_pack_iMul(0,Gz,rhog,grupz);
      mypneb->c_unpack(0,grupx);
      mypneb->c_unpack(0,grupy);
      mypneb->c_unpack(0,grupz);
      mypneb->cr_fft3d(grupx);
      mypneb->cr_fft3d(grupy);
      mypneb->cr_fft3d(grupz);

      /* calculate agrup = |grad nup| */
      mypneb->rr_sqr(grupx,agrup);
      mypneb->rr_addsqr(grupy,agrup);
      mypneb->rr_addsqr(grupz,agrup);
      mypneb->r_sqrt(agrup);


      /* calculate rhodn  */
      mypneb->rr_copy(dn+mypneb->n2ft3d,rhodn);
      mypneb->rr_SMul(scal1,rhodn,rhog);
      mypneb->rc_fft3d(rhog);
      mypneb->c_pack(0,rhog);

      /* calculate   grdn= grad ndn */
      mypneb->tcc_pack_iMul(0,Gx,rhog,grdnx);
      mypneb->tcc_pack_iMul(0,Gy,rhog,grdny);
      mypneb->tcc_pack_iMul(0,Gz,rhog,grdnz);
      mypneb->c_unpack(0,grdnx);
      mypneb->c_unpack(0,grdny);
      mypneb->c_unpack(0,grdnz);
      mypneb->cr_fft3d(grdnx);
      mypneb->cr_fft3d(grdny);
      mypneb->cr_fft3d(grdnz);

      /* calculate agrdn = |grad ndn| */
      mypneb->rr_sqr(grdnx,agrdn);
      mypneb->rr_addsqr(grdny,agrdn);
      mypneb->rr_addsqr(grdnz,agrdn);
      mypneb->r_sqrt(agrdn);

      /* calculate agrall = |grad nup +grad ndn| */
      mypneb->rrr_Sum(grupx,grdnx,grallx);
      mypneb->rrr_Sum(grupy,grdny,grally);
      mypneb->rrr_Sum(grupz,grdnz,grallz);
      mypneb->rr_sqr(grallx,agrall);
      mypneb->rr_addsqr(grally,agrall);
      mypneb->rr_addsqr(grallz,agrall);
      mypneb->r_sqrt(agrall);

      switch (gga) {
         case 10 :
            gen_PBE96_BW_unrestricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 11 :
            gen_BLYP_BW_unrestricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 12 :
            gen_revPBE_BW_unrestricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 13 :
            gen_PBEsol_BW_unrestricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 14 :
            gen_HSE_BW_unrestricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 15 :
            gen_B3LYP_BW_unrestricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
            break;
         case 16 :
            gen_BEEF_BW_unrestricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,0.6001664769,xce,fn,fdn);
            break;
         case 17 :
            gen_BEEF_BW_unrestricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,0.0,xce,fn,fdn);
            break;


         default :
            gen_PBE96_BW_unrestricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);
      }

     /**** calculate df/d|grad nup|* (grad nup)/|grad nup|  ****
      **** calculate df/d|grad ndn|* (grad ndn)/|grad ndn|  ****
      **** calculate df/d|grad n|  * (grad n)/|grad n|  ****/
      mypneb->rr_Divide(agrup,grupx);
      mypneb->rr_Divide(agrup,grupy);
      mypneb->rr_Divide(agrup,grupz);
      mypneb->rr_Divide(agrdn,grdnx);
      mypneb->rr_Divide(agrdn,grdny);
      mypneb->rr_Divide(agrdn,grdnz);
      mypneb->rr_Divide(agrall,grallx);
      mypneb->rr_Divide(agrall,grally);
      mypneb->rr_Divide(agrall,grallz);

      mypneb->rr_Mul(fdnup,grupx);
      mypneb->rr_Mul(fdnup,grupy);
      mypneb->rr_Mul(fdnup,grupz);
      mypneb->rr_Mul(fdndn,grdnx);
      mypneb->rr_Mul(fdndn,grdny);
      mypneb->rr_Mul(fdndn,grdnz);
      mypneb->rr_Mul(fdnall,grallx);
      mypneb->rr_Mul(fdnall,grally);
      mypneb->rr_Mul(fdnall,grallz);

     /**** calculate (df/d|grad nup|* (grad nup)/|grad nup|)  ****
      ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****
      **** calculate (df/d|grad ndn|* (grad ndn)/|grad ndn|)  ****
      ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****/
      mypneb->rr_Sum(grallx,grupx);
      mypneb->rr_Sum(grally,grupy);
      mypneb->rr_Sum(grallz,grupz);
      mypneb->rr_Sum(grallx,grdnx);
      mypneb->rr_Sum(grally,grdny);
      mypneb->rr_Sum(grallz,grdnz);

      mypneb->r_SMul(scal1,grupx);
      mypneb->r_SMul(scal1,grupy);
      mypneb->r_SMul(scal1,grupz);
      mypneb->r_SMul(scal1,grdnx);
      mypneb->r_SMul(scal1,grdny);
      mypneb->r_SMul(scal1,grdnz);

      /* put sums by G-space */
      mypneb->rc_fft3d(grupx);
      mypneb->rc_fft3d(grupy);
      mypneb->rc_fft3d(grupz);
      mypneb->rc_fft3d(grdnx);
      mypneb->rc_fft3d(grdny);
      mypneb->rc_fft3d(grdnz);

      mypneb->c_pack(0,grupx);
      mypneb->c_pack(0,grupy);
      mypneb->c_pack(0,grupz);
      mypneb->c_pack(0,grdnx);
      mypneb->c_pack(0,grdny);
      mypneb->c_pack(0,grdnz);

      /* multiply sums by G vector */
      mypneb->tc_pack_iMul(0,Gx,grupx);
      mypneb->tc_pack_iMul(0,Gy,grupy);
      mypneb->tc_pack_iMul(0,Gz,grupz);
      mypneb->tc_pack_iMul(0,Gx,grdnx);
      mypneb->tc_pack_iMul(0,Gy,grdny);
      mypneb->tc_pack_iMul(0,Gz,grdnz);

      /* addup dot products */
      mypneb->cccc_pack_Sum(0,grupx,grupy,grupz,fdnup);
      mypneb->cccc_pack_Sum(0,grdnx,grdny,grdnz,fdndn);

      /* put back in r-space and subtract from df/dnup,df/dndn */
      mypneb->c_unpack(0,fdnup);
      mypneb->c_unpack(0,fdndn);
      mypneb->cr_fft3d(fdnup);
      mypneb->cr_fft3d(fdndn);
      mypneb->rrr_Minus(fnup,fdnup,xcpup);
      mypneb->rrr_Minus(fndn,fdndn,xcpdn);
   }
}

}

      
#include        "Pneb.hpp"
#include        "pbe96.hpp"
#include        "blyp.hpp"

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
      mypneb->r_nzero(2,agr);

      /* calculate rho tmp1=rho(g) */
      mypneb->rr_copy(dn,rho);
      mypneb->rr_Sum(dn,rho);
      mypneb->rr_SMul(scal1,rho,rhog);
      mypneb->rc_fft3d(rhog);
      mypneb->c_pack(0,rhog);

      /* calculate gr = grad n */
      mypneb->tcc_iMul(0,Gx,rhog,grx);
      mypneb->tcc_iMul(0,Gy,rhog,gry);
      mypneb->tcc_iMul(0,Gz,rhog,grz);
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

      gen_PBE96_BW_restricted(mypneb->n2ft3d,rho,agr,x_parameter,c_parameter,xce,fn,fdn);

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

      mypneb->tc_iMul(0,Gx,grx);
      mypneb->tc_iMul(0,Gy,gry);
      mypneb->tc_iMul(0,Gz,grz);

      mypneb->cccc_Sum(0,grx,gry,grz,fdn);
      mypneb->cr_fft3d(fdn);
      mypneb->rrr_Minus(fn,fdn,xcp);
   } 

  /************************************ 
   ***** unrestricted calculation ***** 
   ************************************/
   else 
   {
   }
}

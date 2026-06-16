
#include "Pneb.hpp"
#include "b3lyp.hpp"
#include "beef_gga.hpp"
#include "blyp.hpp"
#include "hsepbe.hpp"
#include "pbe96.hpp"
#include "pbesol.hpp"
#include "revpbe.hpp"
#include "vdw_DF.hpp"
#include "iofmt.hpp"
#include "units.hpp"

namespace pwdft {

#define dncut 1.0e-30

/********************************
 *				*
 *            v_bwexc           *
 *                              *
 ********************************/
void v_bwexc(const int gga, Pneb *mypneb, vdw_DF *vdw,  const double *dn,
             const double x_parameter, const double c_parameter, double *xcp,
             double *xce, double *rho, double *grx, double *gry, double *grz,
             double *agr, double *fn, double *fdn) 
{
   double *rhog = fn;
   double *Gx = mypneb->Gpackxyz(0, 0);
   double *Gy = mypneb->Gpackxyz(0, 1);
   double *Gz = mypneb->Gpackxyz(0, 2);
   double scal1 = 1.0 / ((double)((mypneb->nx) * (mypneb->ny) * (mypneb->nz)));

   bool has_vdw = false;
   if (vdw != nullptr)
      has_vdw = vdw->exist();
 
   /**********************************
    ***** restricted calculation *****
    **********************************/
   if (mypneb->ispin == 1) 
   {
      // mypneb->r_nzero(2,agr);
      mypneb->r_zero(agr);
     
      /* calculate rho tmp1=rho(g) */
      mypneb->rr_copy(dn, rho);
      mypneb->rr_Sum(dn, rho);
      mypneb->rr_SMul(scal1, rho, rhog);
      mypneb->rc_fft3d(rhog);
      mypneb->c_pack(0, rhog);
     
      /* calculate gr = grad n */
      mypneb->tcc_pack_iMul(0, Gx, rhog, grx);
      mypneb->tcc_pack_iMul(0, Gy, rhog, gry);
      mypneb->tcc_pack_iMul(0, Gz, rhog, grz);
      mypneb->c_unpack(0, grx);
      mypneb->c_unpack(0, gry);
      mypneb->c_unpack(0, grz);
      mypneb->cr_fft3d(grx);
      mypneb->cr_fft3d(gry);
      mypneb->cr_fft3d(grz);
     
      /* calculate agr = |grad n| */
      mypneb->rr_sqr(grx, agr);
      mypneb->rr_addsqr(gry, agr);
      mypneb->rr_addsqr(grz, agr);
      mypneb->r_sqrt(agr);
     
      switch (gga) {
      case 10:
        gen_PBE96_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
        break;
      case 11:
        gen_BLYP_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                               xce, fn, fdn);
        break;
      case 12:
        gen_revPBE_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, xce, fn, fdn);
        break;
      case 13:
        gen_PBEsol_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, xce, fn, fdn);
        break;
      case 14:
        gen_HSE_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                              xce, fn, fdn);
        break;
      case 15:
        gen_B3LYP_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
        break;
      case 16:
        gen_BEEF_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                               0.6001664769, xce, fn, fdn);
        break;
      case 17:
        gen_BEEF_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                               0.0, xce, fn, fdn);
        break;
     
      default:
        gen_PBE96_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);


      }

      // add vdw here
      if (has_vdw)
      {
         vdw->evaluate(mypneb->ispin,dn,agr,xce,fn,fdn);
      }



     
      /* calculate df/d|grad n| *(grad n)/|grad n| */
      mypneb->rr_Divide(agr, grx);
      mypneb->rr_Divide(agr, gry);
      mypneb->rr_Divide(agr, grz);
      mypneb->rr_Mul(fdn, grx);
      mypneb->rr_Mul(fdn, gry);
      mypneb->rr_Mul(fdn, grz);
     
      mypneb->r_SMul(scal1, grx);
      mypneb->r_SMul(scal1, gry);
      mypneb->r_SMul(scal1, grz);
     
      mypneb->rc_fft3d(grx);
      mypneb->rc_fft3d(gry);
      mypneb->rc_fft3d(grz);
     
      mypneb->c_pack(0, grx);
      mypneb->c_pack(0, gry);
      mypneb->c_pack(0, grz);
     
      mypneb->tc_pack_iMul(0, Gx, grx);
      mypneb->tc_pack_iMul(0, Gy, gry);
      mypneb->tc_pack_iMul(0, Gz, grz);
     
      mypneb->cccc_pack_Sum(0, grx, gry, grz, fdn);
      mypneb->c_unpack(0, fdn);
      mypneb->cr_fft3d(fdn);
      mypneb->rrr_Minus(fn, fdn, xcp);
      mypneb->r_zero_ends(xcp);
   }
 
   /************************************
    ***** unrestricted calculation *****
    ************************************/
   else 
   {
      mypneb->r_nzero(3, agr);
     
      double *rhoup = rho;
      double *rhodn = rho + mypneb->n2ft3d;
      double *grupx = grx;
      double *grupy = gry;
      double *grupz = grz;
      double *grdnx = grx + mypneb->n2ft3d;
      double *grdny = gry + mypneb->n2ft3d;
      double *grdnz = grz + mypneb->n2ft3d;
      double *grallx = grx + 2 * mypneb->n2ft3d;
      double *grally = gry + 2 * mypneb->n2ft3d;
      double *grallz = grz + 2 * mypneb->n2ft3d;
      double *agrup = agr;
      double *agrdn = agr + mypneb->n2ft3d;
      double *agrall = agr + 2 * mypneb->n2ft3d;
     
      // double *dnup = dn;
      // double *dndn = dn + mypneb->n2ft3d;
     
      double *xcpup = xcp;
      double *xcpdn = xcp + mypneb->n2ft3d;
     
      double *fnup = fn;
      double *fndn = fn + mypneb->n2ft3d;
      ;
     
      double *fdnup = fdn;
      double *fdndn = fdn + mypneb->n2ft3d;
      ;
      double *fdnall = fdn + 2 * mypneb->n2ft3d;
      ;
     
      /* calculate rhoup  */
      mypneb->rr_copy(dn, rhoup);
      mypneb->rr_SMul(scal1, rhoup, rhog);
      mypneb->rc_fft3d(rhog);
      mypneb->c_pack(0, rhog);
     
      /* calculate   grup= grad nup */
      mypneb->tcc_pack_iMul(0, Gx, rhog, grupx);
      mypneb->tcc_pack_iMul(0, Gy, rhog, grupy);
      mypneb->tcc_pack_iMul(0, Gz, rhog, grupz);
      mypneb->c_unpack(0, grupx);
      mypneb->c_unpack(0, grupy);
      mypneb->c_unpack(0, grupz);
      mypneb->cr_fft3d(grupx);
      mypneb->cr_fft3d(grupy);
      mypneb->cr_fft3d(grupz);
     
      /* calculate agrup = |grad nup| */
      mypneb->rr_sqr(grupx, agrup);
      mypneb->rr_addsqr(grupy, agrup);
      mypneb->rr_addsqr(grupz, agrup);
      mypneb->r_sqrt(agrup);
     
      /* calculate rhodn  */
      mypneb->rr_copy(dn + mypneb->n2ft3d, rhodn);
      mypneb->rr_SMul(scal1, rhodn, rhog);
      mypneb->rc_fft3d(rhog);
      mypneb->c_pack(0, rhog);
     
      /* calculate   grdn= grad ndn */
      mypneb->tcc_pack_iMul(0, Gx, rhog, grdnx);
      mypneb->tcc_pack_iMul(0, Gy, rhog, grdny);
      mypneb->tcc_pack_iMul(0, Gz, rhog, grdnz);
      mypneb->c_unpack(0, grdnx);
      mypneb->c_unpack(0, grdny);
      mypneb->c_unpack(0, grdnz);
      mypneb->cr_fft3d(grdnx);
      mypneb->cr_fft3d(grdny);
      mypneb->cr_fft3d(grdnz);
     
      /* calculate agrdn = |grad ndn| */
      mypneb->rr_sqr(grdnx, agrdn);
      mypneb->rr_addsqr(grdny, agrdn);
      mypneb->rr_addsqr(grdnz, agrdn);
      mypneb->r_sqrt(agrdn);
     
      /* calculate agrall = |grad nup +grad ndn| */
      mypneb->rrr_Sum(grupx, grdnx, grallx);
      mypneb->rrr_Sum(grupy, grdny, grally);
      mypneb->rrr_Sum(grupz, grdnz, grallz);
      mypneb->rr_sqr(grallx, agrall);
      mypneb->rr_addsqr(grally, agrall);
      mypneb->rr_addsqr(grallz, agrall);
      mypneb->r_sqrt(agrall);
     
      switch (gga) {
      case 10:
        gen_PBE96_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                  c_parameter, xce, fn, fdn);
        break;
      case 11:
        gen_BLYP_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, xce, fn, fdn);
        break;
      case 12:
        gen_revPBE_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                   c_parameter, xce, fn, fdn);
        break;
      case 13:
        gen_PBEsol_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                   c_parameter, xce, fn, fdn);
        break;
      case 14:
        gen_HSE_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
        break;
      case 15:
        gen_B3LYP_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                  c_parameter, xce, fn, fdn);
        break;
      case 16:
        gen_BEEF_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, 0.6001664769, xce, fn, fdn);
        break;
      case 17:
        gen_BEEF_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, 0.0, xce, fn, fdn);
        break;
     
      default:
        gen_PBE96_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                  c_parameter, xce, fn, fdn);
      }
     
      // add vdw here
      if (has_vdw)
      {
         std::cout << "vdw HERB" << std::endl;
         vdw->evaluate(mypneb->ispin,dn,agr,xce,fn,fdn);
      }


      /**** calculate df/d|grad nup|* (grad nup)/|grad nup|  ****
       **** calculate df/d|grad ndn|* (grad ndn)/|grad ndn|  ****
       **** calculate df/d|grad n|  * (grad n)/|grad n|  ****/
      mypneb->rr_Divide(agrup, grupx);
      mypneb->rr_Divide(agrup, grupy);
      mypneb->rr_Divide(agrup, grupz);
      mypneb->rr_Divide(agrdn, grdnx);
      mypneb->rr_Divide(agrdn, grdny);
      mypneb->rr_Divide(agrdn, grdnz);
      mypneb->rr_Divide(agrall, grallx);
      mypneb->rr_Divide(agrall, grally);
      mypneb->rr_Divide(agrall, grallz);
     
      mypneb->rr_Mul(fdnup, grupx);
      mypneb->rr_Mul(fdnup, grupy);
      mypneb->rr_Mul(fdnup, grupz);
      mypneb->rr_Mul(fdndn, grdnx);
      mypneb->rr_Mul(fdndn, grdny);
      mypneb->rr_Mul(fdndn, grdnz);
      mypneb->rr_Mul(fdnall, grallx);
      mypneb->rr_Mul(fdnall, grally);
      mypneb->rr_Mul(fdnall, grallz);
     
      /**** calculate (df/d|grad nup|* (grad nup)/|grad nup|)  ****
       ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****
       **** calculate (df/d|grad ndn|* (grad ndn)/|grad ndn|)  ****
       ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****/
      mypneb->rr_Sum(grallx, grupx);
      mypneb->rr_Sum(grally, grupy);
      mypneb->rr_Sum(grallz, grupz);
      mypneb->rr_Sum(grallx, grdnx);
      mypneb->rr_Sum(grally, grdny);
      mypneb->rr_Sum(grallz, grdnz);
     
      mypneb->r_SMul(scal1, grupx);
      mypneb->r_SMul(scal1, grupy);
      mypneb->r_SMul(scal1, grupz);
      mypneb->r_SMul(scal1, grdnx);
      mypneb->r_SMul(scal1, grdny);
      mypneb->r_SMul(scal1, grdnz);
     
      /* put sums by G-space */
      mypneb->rc_fft3d(grupx);
      mypneb->rc_fft3d(grupy);
      mypneb->rc_fft3d(grupz);
      mypneb->rc_fft3d(grdnx);
      mypneb->rc_fft3d(grdny);
      mypneb->rc_fft3d(grdnz);
     
      mypneb->c_pack(0, grupx);
      mypneb->c_pack(0, grupy);
      mypneb->c_pack(0, grupz);
      mypneb->c_pack(0, grdnx);
      mypneb->c_pack(0, grdny);
      mypneb->c_pack(0, grdnz);
     
      /* multiply sums by G vector */
      mypneb->tc_pack_iMul(0, Gx, grupx);
      mypneb->tc_pack_iMul(0, Gy, grupy);
      mypneb->tc_pack_iMul(0, Gz, grupz);
      mypneb->tc_pack_iMul(0, Gx, grdnx);
      mypneb->tc_pack_iMul(0, Gy, grdny);
      mypneb->tc_pack_iMul(0, Gz, grdnz);
     
      /* addup dot products */
      mypneb->cccc_pack_Sum(0, grupx, grupy, grupz, fdnup);
      mypneb->cccc_pack_Sum(0, grdnx, grdny, grdnz, fdndn);
     
      /* put back in r-space and subtract from df/dnup,df/dndn */
      mypneb->c_unpack(0, fdnup);
      mypneb->c_unpack(0, fdndn);
      mypneb->cr_fft3d(fdnup);
      mypneb->cr_fft3d(fdndn);
      mypneb->rrr_Minus(fnup, fdnup, xcpup);
      mypneb->rrr_Minus(fndn, fdndn, xcpdn);
      mypneb->r_zero_ends(xcpup);
      mypneb->r_zero_ends(xcpdn);
   }
}


/********************************
 *				*
 *         v_bwexc_euv          *
 *                              *
 ********************************/
void v_bwexc_euv(const int gga, Pneb *mypneb, vdw_DF *vdw,  const double *dn,
             const double x_parameter, const double c_parameter, 
             double *stress) 
{

   constexpr double pi     = units::PI;
   constexpr double fourpi = 4.0*pi;
   constexpr double scal   = 1.0/(2.0*pi);

   const size_t n2ft3d = mypneb->n2ft3d;
   const double scal1 = 1.0 / ((double)((mypneb->nx) * (mypneb->ny) * (mypneb->nz)));
   const double omega = mypneb->lattice->omega();

   // define hm
   double hm[9];
   for (size_t i=0; i<3; ++i)
   for (size_t j=0; j<3; ++j)
      hm[i+3*j] = scal*mypneb->lattice->unitg(i,j);

   // set G vectors
   const double *Gx = mypneb->Gpackxyz(0,0);
   const double *Gy = mypneb->Gpackxyz(0,1);
   const double *Gz = mypneb->Gpackxyz(0,2);

             
   bool has_vdw = false;
   if (vdw != nullptr)
      has_vdw = vdw->exist();

   std::array<double, 9> W = {0.0}; // Fixed size symmetric tensor
 
   /**********************************
    ***** restricted calculation *****
    **********************************/
   if (mypneb->ispin == 1) 
   {
      // ***** tempory variables needed rho,grx,gry,grz,agr,fn,fdn,rhog,xce *****
      std::vector<double> rho_vec(n2ft3d);
      std::vector<double> grx_vec(n2ft3d);
      std::vector<double> gry_vec(n2ft3d);
      std::vector<double> grz_vec(n2ft3d);
      std::vector<double> agr_vec(n2ft3d);
      std::vector<double> fn_vec(n2ft3d);
      std::vector<double> fdn_vec(2*n2ft3d);
      std::vector<double> rhog_vec(n2ft3d);
      std::vector<double> xce_vec(n2ft3d);
      double *rho = rho_vec.data();
      double *grx = grx_vec.data();
      double *gry = gry_vec.data();
      double *grz = grz_vec.data();
      double *agr = agr_vec.data();
      double *fn  = fn_vec.data();
      double *fdn = fdn_vec.data();
      double *rhog = rhog_vec.data();
      double *xce  = xce_vec.data();

      mypneb->r_zero(agr);
        
      /* calculate rho tmp1=rho(g) */
      mypneb->rr_copy(dn, rho);
      mypneb->rr_Sum(dn, rho);
      mypneb->rr_SMul(scal1, rho, rhog);
      mypneb->rc_fft3d(rhog);   
      mypneb->c_pack(0, rhog);
      
      /* calculate gr = grad n */
      mypneb->tcc_pack_iMul(0, Gx, rhog, grx);
      mypneb->tcc_pack_iMul(0, Gy, rhog, gry);
      mypneb->tcc_pack_iMul(0, Gz, rhog, grz);
      mypneb->c_unpack(0, grx);
      mypneb->c_unpack(0, gry);
      mypneb->c_unpack(0, grz);
      mypneb->cr_fft3d(grx);
      mypneb->cr_fft3d(gry);
      mypneb->cr_fft3d(grz);

      /* calculate agr = |grad n| */
      mypneb->rr_sqr(grx, agr);
      mypneb->rr_addsqr(gry, agr);
      mypneb->rr_addsqr(grz, agr);
      mypneb->r_sqrt(agr);

      switch (gga) {
      case 10:
        gen_PBE96_BW_restricted(n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
        break;
      case 11:
        gen_BLYP_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                               xce, fn, fdn);
        break;
      case 12:
        gen_revPBE_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, xce, fn, fdn);
        break;
      case 13:
        gen_PBEsol_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                 c_parameter, xce, fn, fdn);
        break;
      case 14:
        gen_HSE_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                              xce, fn, fdn);
        break;
      case 15:
        gen_B3LYP_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
        break;
      case 16:
        gen_BEEF_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                               0.6001664769, xce, fn, fdn);
        break;
      case 17:
        gen_BEEF_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter,
                               0.0, xce, fn, fdn);
        break;

      default:
        gen_PBE96_BW_restricted(mypneb->n2ft3d, rho, agr, x_parameter,
                                c_parameter, xce, fn, fdn);
      }

      // add vdw here
      if (has_vdw)
      {
         vdw->evaluate(mypneb->ispin,dn,agr,xce,fn,fdn);
      }

      /* calculate df/d|grad n| *(grad n)/|grad n| */
      mypneb->rr_Divide(agr, grx);
      mypneb->rr_Divide(agr, gry);
      mypneb->rr_Divide(agr, grz);
      mypneb->rr_Mul(fdn, grx);
      mypneb->rr_Mul(fdn, gry);
      mypneb->rr_Mul(fdn, grz);

      mypneb->r_SMul(scal1, grx);
      mypneb->r_SMul(scal1, gry);
      mypneb->r_SMul(scal1, grz);

      mypneb->rc_fft3d(grx);
      mypneb->rc_fft3d(gry);
      mypneb->rc_fft3d(grz);

      mypneb->c_pack(0, grx);
      mypneb->c_pack(0, gry);
      mypneb->c_pack(0, grz);

      mypneb->tc_pack_iMul(0, mypneb->Gpackxyz(0,0), grx);
      mypneb->tc_pack_iMul(0, mypneb->Gpackxyz(0,1), gry);
      mypneb->tc_pack_iMul(0, mypneb->Gpackxyz(0,2), grz);


      //**** W(u,s) = Sum(G) [i*G(u)*dcongj(rhog)*gr(s)] ****
      //****         where gr(1)=grx,gr(2)=gry,gr(3)=grz ****
      for (size_t u=0; u<3; ++u)
      {
         //**** agr = i*G(u)*grx ****
         mypneb->tcc_pack_iMul(0,mypneb->Gpackxyz(0,u),grx,agr);
         double sumx = mypneb->cc_pack_dot(0,rhog,agr);
         W[u] = sumx*omega;

         //**** agr = i*G(u)*gry ****
         mypneb->tcc_pack_iMul(0,mypneb->Gpackxyz(0,u),gry,agr);
         double sumy = mypneb->cc_pack_dot(0,rhog,agr);
         W[u+3] = sumy*omega;

         //**** agr = i*G(u)*grz ****
         mypneb->tcc_pack_iMul(0,mypneb->Gpackxyz(0,u),grz,agr);
         double sumz = mypneb->cc_pack_dot(0,rhog,agr);
         W[u+6] = sumz*omega;
      }

   }
   
   /************************************
    ***** unrestricted calculation *****
    ************************************/
   else 
   {
      // ***** tempory variables needed rho,grx,gry,grz,agr,fn,fdn,rhog,xce *****
      std::vector<double> rho_vec(2*n2ft3d);
      std::vector<double> grxup_vec(n2ft3d);
      std::vector<double> gryup_vec(n2ft3d);
      std::vector<double> grzup_vec(n2ft3d);

      std::vector<double> grxdn_vec(n2ft3d);
      std::vector<double> grydn_vec(n2ft3d);
      std::vector<double> grzdn_vec(n2ft3d);

      std::vector<double> grxall_vec(n2ft3d);
      std::vector<double> gryall_vec(n2ft3d);
      std::vector<double> grzall_vec(n2ft3d);

      std::vector<double> xagr_vec(3*n2ft3d);
      std::vector<double> grad_vec(3*n2ft3d);
      std::vector<double> gtmp_vec(3*n2ft3d);

      std::vector<double> xfn_vec(2*n2ft3d);
      std::vector<double> xfdn_vec(3*n2ft3d);
      std::vector<double> xce_vec(n2ft3d);
      std::vector<double> rhog_vec(2*n2ft3d);

      double *rho = rho_vec.data();
      double *rhoup = rho;
      double *rhodn = rho+n2ft3d;
      double *grupx = grxup_vec.data();
      double *grupy = gryup_vec.data();
      double *grupz = grzup_vec.data();

      double *grdnx = grxdn_vec.data();
      double *grdny = grydn_vec.data();
      double *grdnz = grzdn_vec.data();

      double *grallx = grxall_vec.data();
      double *grally = gryall_vec.data();
      double *grallz = grzall_vec.data();

      double *agr  = xagr_vec.data();
      double *agrup  = agr;
      double *agrdn  = agr+n2ft3d;
      double *agrall = agr+2*n2ft3d;

      double *fn  = xfn_vec.data();
      double *fdn = xfdn_vec.data();

      double *fnup  = fn;
      double *fndn  = fn + n2ft3d;;

      double *fdnup  = fdn;
      double *fdndn  = fdn  + n2ft3d;
      double *fdnall = fdn + 2*n2ft3d;

      double *xce = xce_vec.data();
      double *rhoupg = rhog_vec.data();
      double *rhodng = rhoupg + n2ft3d;


      //***** calculate rhoup  ****
      mypneb->rr_copy(dn, rhoup);
      mypneb->rr_SMul(scal1, rhoup, rhoupg);
      mypneb->rc_fft3d(rhoupg);
      mypneb->c_pack(0, rhoupg);

      //***** calculate   grup= grad nup ****
      mypneb->tcc_pack_iMul(0, Gx, rhoupg, grupx);
      mypneb->tcc_pack_iMul(0, Gy, rhoupg, grupy);
      mypneb->tcc_pack_iMul(0, Gz, rhoupg, grupz);
      mypneb->c_unpack(0, grupx);
      mypneb->c_unpack(0, grupy);
      mypneb->c_unpack(0, grupz);
      mypneb->cr_fft3d(grupx);
      mypneb->cr_fft3d(grupy);
      mypneb->cr_fft3d(grupz);

      //***** calculate agrup = |grad nup| ****
      mypneb->rr_sqr(grupx, agrup);
      mypneb->rr_addsqr(grupy, agrup);
      mypneb->rr_addsqr(grupz, agrup);
      mypneb->r_sqrt(agrup);

      //***** calculate rhodn  ****
      mypneb->rr_copy(dn + mypneb->n2ft3d, rhodn);
      mypneb->rr_SMul(scal1, rhodn, rhodng);
      mypneb->rc_fft3d(rhodng);
      mypneb->c_pack(0, rhodng);

      //***** calculate   grdn= grad ndn ****
      mypneb->tcc_pack_iMul(0, Gx, rhodng, grdnx);
      mypneb->tcc_pack_iMul(0, Gy, rhodng, grdny);
      mypneb->tcc_pack_iMul(0, Gz, rhodng, grdnz);
      mypneb->c_unpack(0, grdnx);
      mypneb->c_unpack(0, grdny);
      mypneb->c_unpack(0, grdnz);
      mypneb->cr_fft3d(grdnx);
      mypneb->cr_fft3d(grdny);
      mypneb->cr_fft3d(grdnz);

      //***** calculate agrdn = |grad ndn| ****
      mypneb->rr_sqr(grdnx, agrdn);
      mypneb->rr_addsqr(grdny, agrdn);
      mypneb->rr_addsqr(grdnz, agrdn);
      mypneb->r_sqrt(agrdn);

      //***** calculate agr = |grad nup +grad ndn| ****
      mypneb->rrr_Sum(grupx, grdnx, grallx);
      mypneb->rrr_Sum(grupy, grdny, grally);
      mypneb->rrr_Sum(grupz, grdnz, grallz);
      mypneb->rr_sqr(grallx, agrall);
      mypneb->rr_addsqr(grally, agrall);
      mypneb->rr_addsqr(grallz, agrall);
      
      switch (gga) {
      case 10:
        gen_PBE96_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 11:
        gen_BLYP_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 12:
        gen_revPBE_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 13:      
        gen_PBEsol_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn); 
        break;
      case 14:
        gen_HSE_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 15:
        gen_B3LYP_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 16:
        gen_BEEF_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter, 0.6001664769, xce, fn, fdn);
        break;
      case 17:
        gen_BEEF_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter, 0.0, xce, fn, fdn);
        break;
      
      default:
        gen_PBE96_BW_unrestricted(mypneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
      }
      
      // add vdw here
      if (has_vdw)
      {
         std::cout << "vdw HERB" << std::endl;
         vdw->evaluate(mypneb->ispin,dn,agr,xce,fn,fdn);
      }

      /**** calculate df/d|grad nup|* (grad nup)/|grad nup|  ****
       **** calculate df/d|grad ndn|* (grad ndn)/|grad ndn|  ****
       **** calculate df/d|grad n|  * (grad n)/|grad n|  ****/
      mypneb->rr_Divide(agrup, grupx);
      mypneb->rr_Divide(agrup, grupy);
      mypneb->rr_Divide(agrup, grupz);
      mypneb->rr_Divide(agrdn, grdnx);
      mypneb->rr_Divide(agrdn, grdny);
      mypneb->rr_Divide(agrdn, grdnz);
      mypneb->rr_Divide(agrall, grallx);
      mypneb->rr_Divide(agrall, grally);
      mypneb->rr_Divide(agrall, grallz);

      mypneb->rr_Mul(fdnup, grupx);
      mypneb->rr_Mul(fdnup, grupy);
      mypneb->rr_Mul(fdnup, grupz);
      mypneb->rr_Mul(fdndn, grdnx);
      mypneb->rr_Mul(fdndn, grdny);
      mypneb->rr_Mul(fdndn, grdnz);
      mypneb->rr_Mul(fdnall, grallx);
      mypneb->rr_Mul(fdnall, grally);
      mypneb->rr_Mul(fdnall, grallz);


      /**** calculate (df/d|grad nup|* (grad nup)/|grad nup|)  ****
       ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****
       **** calculate (df/d|grad ndn|* (grad ndn)/|grad ndn|)  ****
       ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****/
      mypneb->rr_Sum(grallx, grupx);
      mypneb->rr_Sum(grally, grupy);
      mypneb->rr_Sum(grallz, grupz);
      mypneb->rr_Sum(grallx, grdnx);
      mypneb->rr_Sum(grally, grdny);
      mypneb->rr_Sum(grallz, grdnz);

      mypneb->r_SMul(scal1, grupx);
      mypneb->r_SMul(scal1, grupy);
      mypneb->r_SMul(scal1, grupz);
      mypneb->r_SMul(scal1, grdnx);
      mypneb->r_SMul(scal1, grdny);
      mypneb->r_SMul(scal1, grdnz);

      /* put sums by G-space */
      mypneb->rc_fft3d(grupx);
      mypneb->rc_fft3d(grupy);
      mypneb->rc_fft3d(grupz);
      mypneb->rc_fft3d(grdnx);
      mypneb->rc_fft3d(grdny);
      mypneb->rc_fft3d(grdnz);

      mypneb->c_pack(0, grupx);
      mypneb->c_pack(0, grupy);
      mypneb->c_pack(0, grupz);
      mypneb->c_pack(0, grdnx);
      mypneb->c_pack(0, grdny);
      mypneb->c_pack(0, grdnz);

      /* multiply sums by G vector */
      mypneb->tc_pack_iMul(0, Gx, grupx);
      mypneb->tc_pack_iMul(0, Gy, grupy);
      mypneb->tc_pack_iMul(0, Gz, grupz);
      mypneb->tc_pack_iMul(0, Gx, grdnx);
      mypneb->tc_pack_iMul(0, Gy, grdny);
      mypneb->tc_pack_iMul(0, Gz, grdnz);

      //**** W(u,s) = Sum(G) [i*G(u)*dcongj(rhoup)*grup(s)] ****
      //****         where grup(1)=grupx,grup(2)=grupy,grup(3)=grupz ****
      for (size_t u=0; u<3; ++u)
      {
         //**** agr = i*G(u)*grx ****
         mypneb->tcc_pack_iMul(0,mypneb->Gpackxyz(0,u),grupx,agrup);
         double sumx = mypneb->cc_pack_dot(0,rhoupg,agrup);
         W[u] = sumx*omega;

         //**** agr = i*G(u)*gry ****
         mypneb->tcc_pack_iMul(0,mypneb->Gpackxyz(0,u),grupy,agrup);
         double sumy = mypneb->cc_pack_dot(0,rhoupg,agrup);
         W[u+3] = sumy*omega;

         //**** agr = i*G(u)*grz ****
         mypneb->tcc_pack_iMul(0,mypneb->Gpackxyz(0,u),grupz,agrup);
         double sumz = mypneb->cc_pack_dot(0,rhoupg,agrup);
         W[u+6] = sumz*omega;
      }

      //**** W(u,s) = Sum(G) [i*G(u)*dcongj(rhodn)*grup(s)] ****
      //****         where grdn(1)=grdnx,grdn(2)=grdny,grdn(3)=grdnz ****
      for (size_t u=0; u<3; ++u)
      {
         //**** agr = i*G(u)*grx ****
         mypneb->tcc_pack_iMul(0,mypneb->Gpackxyz(0,u),grdnx,agrdn);
         double sumx = mypneb->cc_pack_dot(0,rhodng,agrdn);
         W[u] += sumx*omega;

         //**** agr = i*G(u)*gry ****
         mypneb->tcc_pack_iMul(0,mypneb->Gpackxyz(0,u),grdny,agrdn);
         double sumy = mypneb->cc_pack_dot(0,rhodng,agrdn);
         W[u+3] += sumy*omega;

         //**** agr = i*G(u)*grz ****
         mypneb->tcc_pack_iMul(0,mypneb->Gpackxyz(0,u),grdnz,agrdn);
         double sumz = mypneb->cc_pack_dot(0,rhodng,agrdn);
         W[u+6] += sumz*omega;
      }


      


   }

   for (size_t v=0; v<3; ++v)
   for (size_t u=0; u<3; ++u)
     stress[u+3*v] += W[u+3*v];

}

} // namespace pwdft

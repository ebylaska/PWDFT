
#include "Cneb.hpp"
#include "b3lyp.hpp"
#include "beef_gga.hpp"
#include "blyp.hpp"
#include "hsepbe.hpp"
#include "pbe96.hpp"
#include "pbesol.hpp"
#include "revpbe.hpp"
#include "cvdw_DF.hpp"
#include "units.hpp"

namespace pwdft {

#define dncut 1.0e-30

/********************************
 *				                    *
 *            v_cwexc           *
 *                              *
 ********************************/
// rho,agr,fn,fdn - real
// grx,gry,grz - complex
void v_cwexc(const int gga, Cneb *mycneb,  cvdw_DF *vdw, const double *dn,
             const double x_parameter, const double c_parameter, double *xcp,
             double *xce, double *rho, double *grx, double *gry, double *grz,
             double *agr, double *fn, double *fdn) 
{
   double *rhog = fn; // complex
   double *Gx = mycneb->Gpackxyz(0, 0);
   double *Gy = mycneb->Gpackxyz(0, 1);
   double *Gz = mycneb->Gpackxyz(0, 2);
   double scal1 = 1.0 / ((double)((mycneb->nx) * (mycneb->ny) * (mycneb->nz)));

   bool has_vdw = false;
   if (vdw != nullptr)
      has_vdw = vdw->exist();
 
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

      // add vdw here
      if (has_vdw)
      {
         vdw->evaluate(mycneb->ispin,dn,agr,xce,fn,fdn);
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

      // add vdw here
      if (has_vdw)
      {
         vdw->evaluate(mycneb->ispin,dn,agr,xce,fn,fdn);
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


/********************************
 *				                    *
 *         v_cwexc_euv          *
 *                              *
 ********************************/
void v_cwexc_euv(const int gga, Cneb *mycneb, cvdw_DF *vdw,  const double *dn,
             const double x_parameter, const double c_parameter,
             double *stress)
{

   constexpr double pi     = units::PI;
   constexpr double fourpi = 4.0*pi;
   constexpr double scal   = 1.0/(2.0*pi);

   const size_t n2ft3d = mycneb->n2ft3d;
   const double scal1 = 1.0 / ((double)((mycneb->nx) * (mycneb->ny) * (mycneb->nz)));
   const double omega = mycneb->lattice->omega();

   // define hm
   double hm[9];
   for (size_t i=0; i<3; ++i)
   for (size_t j=0; j<3; ++j)
      hm[i+3*j] = scal*mycneb->lattice->unitg(i,j);

   // set G vectors
   const double *Gx = mycneb->Gpackxyz(0,0);
   const double *Gy = mycneb->Gpackxyz(0,1);
   const double *Gz = mycneb->Gpackxyz(0,2);


   bool has_vdw = false;
   if (vdw != nullptr)
      has_vdw = vdw->exist();

   std::array<double, 9> W = {0.0}; // Fixed size symmetric tensor

   /**********************************
    ***** restricted calculation *****
    **********************************/
   if (mycneb->ispin == 1)
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

      mycneb->r_zero(agr);

      /*ccalculate rho tmp1=rho(g) */
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
        gen_PBE96_BW_restricted(n2ft3d, rho, agr, x_parameter,
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

      // add vdw here
      if (has_vdw)
      {
         vdw->evaluate(mycneb->ispin,dn,agr,xce,fn,fdn);
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


      //mypneb->tc_pack_iMul(0, mypneb->Gpackxyz(0,0), grx);
      //mypneb->tc_pack_iMul(0, mypneb->Gpackxyz(0,1), gry);
      //mypneb->tc_pack_iMul(0, mypneb->Gpackxyz(0,2), grz);

      //**** W(u,s) = Sum(G) [i*G(u)*dcongj(rhog)*gr(s)] ****
      //****         where gr(1)=grx,gr(2)=gry,gr(3)=grz ****
      for (size_t u=0; u<3; ++u)
      {
         //**** agr = i*G(u)*grx ****
         mycneb->tcc_pack_iMul(0,mycneb->Gpackxyz(0,u),grx,agr);
         double sumx = mycneb->cc_pack_dot(0,rhog,agr);
         W[u] = sumx*omega;

         //**** agr = i*G(u)*gry ****
         mycneb->tcc_pack_iMul(0,mycneb->Gpackxyz(0,u),gry,agr);
         double sumy = mycneb->cc_pack_dot(0,rhog,agr);
         W[u+3] = sumy*omega;

         //**** agr = i*G(u)*grz ****
         mycneb->tcc_pack_iMul(0,mycneb->Gpackxyz(0,u),grz,agr);
         double sumz = mycneb->cc_pack_dot(0,rhog,agr);
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
      mycneb->rr_copy(dn, rhoup);
      mycneb->rr_SMul(scal1, rhoup, rhoupg);
      mycneb->rc_fft3d(rhoupg);
      mycneb->c_pack(0, rhoupg);

      //***** calculate   grup= grad nup ****
      mycneb->tcc_pack_iMul(0, Gx, rhoupg, grupx);
      mycneb->tcc_pack_iMul(0, Gy, rhoupg, grupy);
      mycneb->tcc_pack_iMul(0, Gz, rhoupg, grupz);
      mycneb->c_unpack(0, grupx);
      mycneb->c_unpack(0, grupy);
      mycneb->c_unpack(0, grupz);
      mycneb->cr_fft3d(grupx);
      mycneb->cr_fft3d(grupy);
      mycneb->cr_fft3d(grupz);

      //***** calculate agrup = |grad nup| ****
      mycneb->rr_sqr(grupx, agrup);
      mycneb->rr_addsqr(grupy, agrup);
      mycneb->rr_addsqr(grupz, agrup);
      mycneb->r_sqrt(agrup);

      //***** calculate rhodn  ****
      mycneb->rr_copy(dn + mycneb->n2ft3d, rhodn);
      mycneb->rr_SMul(scal1, rhodn, rhodng);
      mycneb->rc_fft3d(rhodng);
      mycneb->c_pack(0, rhodng);

      //***** calculate   grdn= grad ndn ****
      mycneb->tcc_pack_iMul(0, Gx, rhodng, grdnx);
      mycneb->tcc_pack_iMul(0, Gy, rhodng, grdny);
      mycneb->tcc_pack_iMul(0, Gz, rhodng, grdnz);
      mycneb->c_unpack(0, grdnx);
      mycneb->c_unpack(0, grdny);
      mycneb->c_unpack(0, grdnz);
      mycneb->cr_fft3d(grdnx);
      mycneb->cr_fft3d(grdny);
      mycneb->cr_fft3d(grdnz);

      //***** calculate agrdn = |grad ndn| ****
      mycneb->rr_sqr(grdnx, agrdn);
      mycneb->rr_addsqr(grdny, agrdn);
      mycneb->rr_addsqr(grdnz, agrdn);
      mycneb->r_sqrt(agrdn);

      //***** calculate agr = |grad nup +grad ndn| ****
      mycneb->rrr_Sum(grupx, grdnx, grallx);
      mycneb->rrr_Sum(grupy, grdny, grally);
      mycneb->rrr_Sum(grupz, grdnz, grallz);
      mycneb->rr_sqr(grallx, agrall);
      mycneb->rr_addsqr(grally, agrall);
      mycneb->rr_addsqr(grallz, agrall);

      switch (gga) {
      case 10:
        gen_PBE96_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 11:
        gen_BLYP_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 12:
        gen_revPBE_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 13:
        gen_PBEsol_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 14:
        gen_HSE_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 15:
        gen_B3LYP_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
        break;
      case 16:
        gen_BEEF_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter, 0.6001664769, xce, fn, fdn);
        break;
      case 17:
        gen_BEEF_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter, 0.0, xce, fn, fdn);
        break;

      default:
        gen_PBE96_BW_unrestricted(mycneb->n2ft3d, rho, agr, x_parameter, c_parameter, xce, fn, fdn);
      }

      // add vdw here
      if (has_vdw)
      {
         std::cout << "vdw HERB" << std::endl;
         vdw->evaluate(mycneb->ispin,dn,agr,xce,fn,fdn);
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
      mycneb->tc_pack_iMul(0, Gx, grupx);
      mycneb->tc_pack_iMul(0, Gy, grupy);
      mycneb->tc_pack_iMul(0, Gz, grupz);
      mycneb->tc_pack_iMul(0, Gx, grdnx);
      mycneb->tc_pack_iMul(0, Gy, grdny);
      mycneb->tc_pack_iMul(0, Gz, grdnz);

      //**** W(u,s) = Sum(G) [i*G(u)*dcongj(rhoup)*grup(s)] ****
      //****         where grup(1)=grupx,grup(2)=grupy,grup(3)=grupz ****
      for (size_t u=0; u<3; ++u)
      {
         //**** agr = i*G(u)*grx ****
         mycneb->tcc_pack_iMul(0,mycneb->Gpackxyz(0,u),grupx,agrup);
         double sumx = mycneb->cc_pack_dot(0,rhoupg,agrup);
         W[u] = sumx*omega;

         //**** agr = i*G(u)*gry ****
         mycneb->tcc_pack_iMul(0,mycneb->Gpackxyz(0,u),grupy,agrup);
         double sumy = mycneb->cc_pack_dot(0,rhoupg,agrup);
         W[u+3] = sumy*omega;

         //**** agr = i*G(u)*grz ****
         mycneb->tcc_pack_iMul(0,mycneb->Gpackxyz(0,u),grupz,agrup);
         double sumz = mycneb->cc_pack_dot(0,rhoupg,agrup);
         W[u+6] = sumz*omega;
      }

      //**** W(u,s) = Sum(G) [i*G(u)*dcongj(rhodn)*grup(s)] ****
      //****         where grdn(1)=grdnx,grdn(2)=grdny,grdn(3)=grdnz ****
      for (size_t u=0; u<3; ++u)
      {
         //**** agr = i*G(u)*grx ****
         mycneb->tcc_pack_iMul(0,mycneb->Gpackxyz(0,u),grdnx,agrdn);
         double sumx = mycneb->cc_pack_dot(0,rhodng,agrdn);
         W[u] += sumx*omega;

         //**** agr = i*G(u)*gry ****
         mycneb->tcc_pack_iMul(0,mycneb->Gpackxyz(0,u),grdny,agrdn);
         double sumy = mycneb->cc_pack_dot(0,rhodng,agrdn);
         W[u+3] += sumy*omega;

         //**** agr = i*G(u)*grz ****
         mycneb->tcc_pack_iMul(0,mycneb->Gpackxyz(0,u),grdnz,agrdn);
         double sumz = mycneb->cc_pack_dot(0,rhodng,agrdn);
         W[u+6] += sumz*omega;
      }





   }


   for (size_t v=0; v<3; ++v)
   {
      for (size_t u=0; u<3; ++u)
      {
         stress[u+3*v] = 0.0;
         for (size_t s=0; s<3; ++s)
            stress[u+3*v] += W[u+3*s]*hm[s+3*v];
      }
   }


}

} // namespace pwdft

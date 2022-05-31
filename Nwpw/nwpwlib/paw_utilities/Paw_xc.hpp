#ifndef _PAW_XC_HPP_
#define _PAW_XC_HPP_

#include        "Ion.hpp"
#include        "Ewald.hpp"
#include	"Pneb.hpp"
#include        "Control2.hpp"
#include	"util_gaunt.hpp"

namespace pwdft {


class Paw_xc {

   Ion      *myion;
   Pneb     *mypneb;

   int mult_l_max,gga;
   int paw_xc_work_size,n1dgrid_max,nbasis_max,nprj_max;
   int paw_xc_nphi, paw_xc_ntheta, paw_xc_pot_size;

   int *nindx_rholm, *shift_rholm;
   int *bi_rholm, *bj_rholm, *iprj_rholm, *jprj_rholm, *lm_rholm;

   double *coeff_rholm, *coeff_rholm2, *coeff_rholm3;

   double *paw_rho2_ae, *paw_rho2_ps, *paw_rho2_ae_prime, *paw_rho2_ps_prime;
   double *rho_ae, *rho_ps, *rho_ae_prime, *rho_ps_prime, *agr_ae, *agr_ps;
   double *vxc_ae, *vxc_ps, *fdn_ae, *fdn_ps, *exc_ae, *exc_ps, *xc_temp, *xc_temp2;

   double *paw_xc_tlm, *paw_xc_dtlm_theta, *paw_xc_dtlm_phi;
   double *paw_xc_angle_phi, *paw_xc_cos_theta, *paw_xc_w_phi, *paw_xc_w_theta,*paw_xc_e;
   double *paw_vxc_ae, *paw_vxc_ps, *paw_dvxc_ae, *paw_dvxc_ps, *paw_xc_pot, *paw_xc_matr, *paw_xc_dmatr;



public:

   /* Constructors */
   Paw_xc(Ion *,Pneb *,Control2&,
          const int *, const int *, const int *, const int *, const int *,
          int **, int **, int **);


   /* destructor */
   ~Paw_xc() {

      if (gga>=10)
      {
         delete [] paw_rho2_ae_prime;
         delete [] paw_rho2_ps_prime;
         delete [] rho_ae_prime;
         delete [] rho_ps_prime; 
         delete [] agr_ae;
         delete [] agr_ps;
         delete [] fdn_ae;
         delete [] fdn_ps;
         delete [] paw_xc_dtlm_theta;
         delete [] paw_xc_dtlm_phi;
         delete [] paw_xc_dmatr;
         delete [] paw_dvxc_ae; 
         delete [] paw_dvxc_ps;
      }

      delete [] paw_rho2_ae;
      delete [] paw_rho2_ps;
      delete [] rho_ae;
      delete [] rho_ps;

      delete [] vxc_ae;
      delete [] vxc_ps;
      delete [] exc_ae;
      delete [] exc_ps;
      delete [] xc_temp;
      delete [] xc_temp2;

      delete [] paw_xc_tlm;
      delete [] paw_xc_angle_phi;
      delete [] paw_xc_cos_theta;
      delete [] paw_xc_w_phi;
      delete [] paw_xc_w_theta;
      delete [] paw_xc_e;
      delete [] paw_vxc_ae;
      delete [] paw_vxc_ps; 
      delete [] paw_xc_pot; 
      delete [] paw_xc_matr;

      delete [] nindx_rholm;
      delete [] shift_rholm;
      delete [] bi_rholm;
      delete [] bj_rholm;
      delete [] iprj_rholm;
      delete [] jprj_rholm;
      delete [] lm_rholm;
      delete [] coeff_rholm;
      delete [] coeff_rholm2;
      delete [] coeff_rholm3;

   }

   /* sets the gaussian integrals */
   //void set(const bool);

   /* sets the ewald */
   //void set_ewald(Ewald *myewald) { myewald = myewald0; }

};

}

#endif


#pragma once

//

#include        "gdevice.hpp"
#include	"Control2.hpp"
#include	"Ion.hpp"
#include	"Pneb.hpp"
#include	"Strfac.hpp"
#include	"Paw_compcharge.hpp"
#include	"Paw_xc.hpp"
#include	"nwpw_apc.hpp"

namespace pwdft {

class	Pseudopotential {

   int nprj_max;

   double **Gijl;
   double **ncore_atom;
   double **vl;
   double **vlpaw;
   double **vnl;
   double *amass;

   //char **atomsym;
   Pneb   *mypneb;
   Ion    *myion;
   Strfac *mystrfac;

public:
   nwpw_apc       *myapc;
   Paw_compcharge *mypaw_compcharge;
   Paw_xc         *mypaw_xc;

   int psp_version = 3;

   bool *semicore;
   int npsp;
   int *nprj,*lmax,*lmmax,*locp,*nmax,*psp_type;

   int **n_projector,**l_projector,**m_projector,**b_projector;
   double **rc;
   double *zv,*rlocal,*rcore,*ncore_sum;
   double *semicore_density;
   char **comment;

   // paw variables
   bool pawexist=false;
   double **hartree_matrix,**comp_charge_matrix,**comp_pot_matrix;
   double **rgrid,**eig,**phi_ae,**dphi_ae, **phi_ps,**dphi_ps;
   double **core_ae,**core_ps,**core_ae_prime,**core_ps_prime;

   double *log_amesh,*r1,*rmax,*sigma,*zion,*core_kin,*core_ion;
   int    *n1dgrid,*n1dbasis,**nae,**nps,**lps,*icut;


   /* Constructors */
   Pseudopotential(Ion *, Pneb *, Strfac *, Control2&, std::ostream&);

   /* destructor */
   ~Pseudopotential() {
      for (int ia=0; ia<npsp; ++ia)
      {
         delete [] vl[ia];
         delete [] vlpaw[ia];
         delete [] comment[ia];
         delete [] rc[ia];
         if (nprj[ia]>0)
         {
            delete [] n_projector[ia];
            delete [] l_projector[ia];
            delete [] m_projector[ia];
            delete [] b_projector[ia];
            delete [] Gijl[ia];
            delete [] vnl[ia];
         }
         if (semicore[ia])
            delete [] ncore_atom[ia];

         if (psp_type[ia]==4)
         {
            delete [] nae[ia];
            delete [] nps[ia];
            delete [] lps[ia];
            delete [] eig[ia];
            delete [] phi_ae[ia];
            delete [] dphi_ae[ia];
            delete [] phi_ps[ia];
            delete [] dphi_ps[ia];
            delete [] core_ae[ia];
            delete [] core_ps[ia];
            delete [] core_ae_prime[ia];
            delete [] core_ps_prime[ia];
            delete [] rgrid[ia];

            delete [] hartree_matrix[ia];
            delete [] comp_charge_matrix[ia];
            delete [] comp_pot_matrix[ia];
         }

      }
      delete [] vl;
      delete [] vnl;
      delete [] ncore_atom;
      delete [] comment;
      delete [] nprj;
      delete [] lmax;
      delete [] locp;
      delete [] nmax;
      delete [] psp_type;
      delete [] zv;
      delete [] n_projector;
      delete [] l_projector;
      delete [] m_projector;
      delete [] b_projector;
      delete [] rlocal;
      delete [] semicore;
      delete [] rcore;
      delete [] ncore_sum;
      delete [] rc;

      // *** paw data ***
      delete [] hartree_matrix;
      delete [] comp_charge_matrix;
      delete [] comp_pot_matrix;
      delete [] vlpaw;

      delete [] log_amesh; 
      delete [] r1;
      delete [] rmax;
      delete [] sigma;
      delete [] zion;
      delete [] core_kin;
      delete [] core_ion;
      delete [] n1dgrid;
      delete [] n1dbasis;
      delete [] nae;
      delete [] nps;
      delete [] lps;
      delete [] icut;
      delete [] eig;
      delete [] phi_ae;
      delete [] dphi_ae;
      delete [] phi_ps;
      delete [] dphi_ps;
      delete [] core_ae;
      delete [] core_ps;
      delete [] core_ae_prime;
      delete [] core_ps_prime;
      delete [] rgrid;

      if (pawexist) 
      {
         delete mypaw_compcharge;
         delete mypaw_xc;
      }

      delete myapc;
      mypneb->r_dealloc(semicore_density);
    }

    bool has_semicore() { return semicore[npsp]; }
    double ncore(const int ia) {return ncore_sum[ia];}
    void semicore_density_update();
    void semicore_xc_fion(double *, double *);

    void v_nonlocal(double *, double *);
    void v_nonlocal_fion(double *, double *, const bool, double *);
    void f_nonlocal_fion(double *, double *);

    void v_local(double *, const bool, double *, double *);
    void f_local(double *, double *);

    void v_lr_local(double *);
    void grad_v_lr_local(const double *, double *);

    double e_nonlocal(double *);

    double sphere_radius(const int ia) { return rgrid[ia][icut[ia]-1]; }

    std::string print_pspall();

    char spdf_name(const int l) {
       char name = '?';
       if (l==0) name = 's';
       if (l==1) name = 'p';
       if (l==2) name = 'd';
       if (l==3) name = 'f';
       if (l==4) name = 'g';
       if (l==5) name = 'h';
       if (l==6) name = 'i';
       if (l==7) name = 'j';
       if (l==8) name = 'k';
       if (l==9) name = 'l';
       return name;
    }

};
 
}

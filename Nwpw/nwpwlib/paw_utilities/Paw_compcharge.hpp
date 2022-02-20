#ifndef _PAW_COMPCHARGE_HPP_
#define _PAW_COMPCHARGE_HPP_

#include        "Ion.hpp"
#include        "Ewald.hpp"
#include	"Pneb.hpp"
#include        "Control2.hpp"
#include	"util_gaunt.hpp"

namespace pwdft {
using namespace pwdft;

class Paw_compcharge {

   Ion      *myion;
   Pneb     *mypneb;

   Ewald    *myewald;

   //**** common block for nwpw_compcharge data  ****
   bool  isgamma,use_grid_cmp;
   int mult_l_max,lm_size_max; 
   int nion_paw,nkatm_paw,*katm_paw;
   int *katm_pawtoion,*ion_pawtoion;
   int *katm_iontopaw,*ion_iontopaw;
   int *mult_l,*lm_size;

   double *sigma_paw,sigma_smooth;

   double *vk_smooth;
   double *gk_smooth;     // gk_smooth(k)  = 4*pi * Exp[-k*k*sigma_smooth**2 / 4] 
   double *gk;           // gk(k,1:nkatm_paw) = 4*pi * Exp[-k*k*sigma(ia)**2 / 4]
   double *glm;          // glm(k, lm=1:(max_mult_l+1)**2) =  Tlm(k) * |k|**l / (2*l+1)!!
                          //    - Note that (-i)**l factor will be assumed in the calculation.
                          //      Also note that  the Tlm and gaussian terms are rearranged to save space
                          //      It is more natural to define g_lm, gk and gk_smooth
                          //           g_lm      =  4*pi*Tlm(k)
                          //           gk        = |k|**l * Exp[-k*k*sigma(ia)**2 / 4] / (2*l+1)!!
                          //           gk_smooth = |k|**l * Exp[-k*k*sigma_smooth**2 / 4] / (2*l+1)!!
                          //       but this requires more space

   double *Qlm;          // Qlm(lm=1:(mult_l+1)**2,2,1:nion_paw) = compensation charge coefficients
   double *Qlmx;         // Qlmx(lm=1:(mult_l+1)**2,2,1:nion_paw) = compensation charge coefficients
   double *Qlmy;         // Qlmy(lm=1:(mult_l+1)**2,2,1:nion_paw) = compensation charge coefficients
   double *Qlmz;         // Qlmz(lm=1:(mult_l+1)**2,2,1:nion_paw) = compensation charge coefficients
   double *dEmult_Qlm;   // dEmult_Qlm(lm=1:(mult_l+1)**2,2,1:nion_paw) = dEmult/dQlm
   double *dElocal_Qlm;  // dElocal_Qlm(lm=1:(mult_l+1)**2,2,1:nion_paw) = dElocal/dQlm
   double *dE_Qlm;       // dE_Qlm(lm=1:(mult_l+1)**2,2,1:nion_paw) = dE/dQlm where E is the total energy of the system

   int *nindx_Tndiff;
   int *shift_Tndiff;
   int *lm_Tndiff;
   int *iprj_Tndiff;
   int *jprj_Tndiff;
   double *coeff_Tndiff;

   //**** data for atomic hartree indexing ****
   int *nindx_hartree;
   int *shift_hartree;
   int *iprj_hartree;
   int *jprj_hartree;
   int *iprj1_hartree;
   int *jprj1_hartree;
   double *coeff_hartree;


public:

   /* Constructors */
   Paw_compcharge(Ion *, Pneb *, Control2&,
                  const int *, const int *, const int *, const int *,
                  const double *,
                  const int, int **, int **, int **,
                  double **, double **);

   /* destructor */
   ~Paw_compcharge() {

      util_gaunt_end();

      delete [] vk_smooth;
      delete [] gk_smooth;

      delete [] glm;

      delete [] gk;
      delete [] Qlm;
      delete [] Qlmx;
      delete [] Qlmy;
      delete [] Qlmz;
      delete [] dEmult_Qlm;
      delete [] dElocal_Qlm;
      delete [] dE_Qlm;
   }

   /* sets the gaussian integrals */
   //void set(const bool);

   /* sets the ewald */
   //void set_ewald(Ewald *myewald) { myewald = myewald0; }

};

}

#endif


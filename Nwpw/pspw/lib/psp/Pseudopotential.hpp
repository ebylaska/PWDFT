#ifndef	_PSEUDOPOTENTIAL_H_
#define _PSEUDOPOTENTIAL_H_

using namespace std;

#pragma once

#include        "gdevice.hpp"
#include	"Control2.hpp"
#include	"Ion.hpp"
#include	"Pneb.hpp"
#include	"Strfac.hpp"

class	Pseudopotential {

   int nprj_max;

   double **Gijl;
   double **ncore_atom;
   double **vl;
   double **vnl;
   double *rlocal;
   double *amass;

   //char **atomsym;
   Pneb   *mypneb;
   Ion    *myion;
   Strfac *mystrfac;

public:
   bool *semicore;
   int npsp;
   int *nprj,*lmax,*lmmax,*locp,*nmax,*psp_type;

   int **n_projector,**l_projector,**m_projector,**b_projector;
   double **rc;
   double *zv,*rcore,*ncore_sum;
   double *semicore_density;
   char **comment;

   /* Constructors */
   Pseudopotential(Ion *, Pneb *, Strfac *, Control2&);

   /* destructor */
   ~Pseudopotential() {
      for (int ia=0; ia<npsp; ++ia)
      {
         delete [] vl[ia];
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
      delete [] semicore;
      delete [] rcore;
      delete [] ncore_sum;
      delete [] rc;

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

    double e_nonlocal(double *);

};

#endif

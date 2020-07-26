#ifndef	_PSP1D_HAMANN_H_
#define _PSP1D_HAMANN_H_

#include        "Parallel.hpp"

class	Psp1d_Hamann {

public:
   int version,nrho,nmax,lmax0,lmax,locp,psp_type;
   int nprj,n_extra,n_expansion[10];
   double drho,rlocal,amass,zv;

   double rc[10];
   char atom[2];
   char comment[80];

   double *rho;
   double *vp;
   double *wp;
   double *vnlnrm;

   double *up;
   double *r3_matrix;

   bool   semicore;
   int    isemicore;
   double rcore;
   double *rho_sc_r;


   /* Constructors */
   Psp1d_Hamann(Parallel *, const char *);

   /* destructor */
   ~Psp1d_Hamann() {
      delete [] rho;
      delete [] vp;
      delete [] wp;
      delete [] vnlnrm;
      if (psp_type==9)
      {
         delete [] up;
         delete [] r3_matrix;
      }

      if (semicore) delete [] rho_sc_r;
    }

   /* G integration routines */
   void vpp_generate_ray(Parallel *, int, double *, double *, double *, double *);
};

#endif

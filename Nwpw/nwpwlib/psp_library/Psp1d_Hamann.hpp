#ifndef	_PSP1D_HAMANN_H_
#define _PSP1D_HAMANN_H_

#include        "Parallel.hpp"
#include        "PGrid.hpp"

class	Psp1d_Hamann {

public:
   int version,nrho,nmax,lmax0,lmax,locp,psp_type;
   int nprj,n_extra,n_expansion[10];
   float drho,rlocal,amass,zv;

   float rc[10];
   char atom[2];
   char comment[80];

   float *rho;
   float *vp;
   float *wp;
   float *vnlnrm;
   int    *n_prj, *l_prj, *m_prj, *b_prj;

   float *up;
   float *r3_matrix;

   bool   semicore;
   int    isemicore;
   float rcore;
   float *rho_sc_r;


   /* Constructors */
   Psp1d_Hamann(Parallel *, const char *);

   /* destructor */
   ~Psp1d_Hamann() {
      delete [] rho;
      delete [] vp;
      delete [] wp;
      delete [] vnlnrm;
      if (nprj>0)
      {
         delete [] n_prj;
         delete [] l_prj;
         delete [] m_prj;
         delete [] b_prj;
      }
      if (psp_type==9)
      {
         delete [] up;
         delete [] r3_matrix;
      }

      if (semicore) delete [] rho_sc_r;
    }

   /* G integration routines */
   void vpp_generate_ray(Parallel *, int, float *, float *, float *, float *);

   void vpp_generate_spline(PGrid *, int, float *, float *, float *, float *, float *, float *, float *);
};

#endif

#ifndef	_PSP1D_HAMANN_H_
#define _PSP1D_HAMANN_H_

#include        "Parallel.hpp"

class	Psp1d_Hamann {

public:
   int version,nrho,nmax,lmax0,lmax,locp,ihasae;
   double drho,rlocal,amass,zv;

   double rc[10];
   char atom[2];
   char comment[80];

   double *rho;
   double *vp;
   double *wp;

   double *up;
   double *r3_matrix;

   int    isemicore;
   double rcore;
   double *semicore;


   /* Constructors */
   Psp1d_Hamann(Parallel *, const char *);

   /* destructor */
   ~Psp1d_Hamann() {
      delete [] rho;
      delete [] vp;
      delete [] wp;
      if (ihasae>0)
      {
         delete [] up;
         delete [] r3_matrix;
      }
      if (isemicore>0) delete [] semicore;
    }
};

#endif

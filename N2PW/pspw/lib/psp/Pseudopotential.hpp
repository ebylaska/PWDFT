#ifndef	_PSEUDOPOTENTIAL_H_
#define _PSEUDOPOTENTIAL_H_




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
   int npsp;
   int *nprj,*lmax,*lmmax,*locp,*nmax,*psp_type,*semicore;
   int **n_projector,**l_projector,**m_projector,**b_projector;
   double **rc;
   double *zv,*rcore;
   char **comment;

   /* Constructors */
   Pseudopotential(Ion *, Pneb *, Strfac *);

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
      delete [] semicore;
      delete [] rcore;
   //   delete [] vl;
   //   delete [] vnl;
   //   delete [] rlocal;
      delete [] rc;
    }

    double ncore(const int ia) {return 0.0;}

    void v_nonlocal(double *, double *);
    void v_local(double *);

};

#endif

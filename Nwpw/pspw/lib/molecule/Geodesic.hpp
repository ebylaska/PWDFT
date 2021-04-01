#ifndef	_GEODESIC_HPP_
#define _GEODESIC_HPP_

#include	<cmath>
#include	"Pneb.hpp"
#include	"Electron.hpp"
#include	"Molecule.hpp"

class	Geodesic {

   int minimizer;
   Pneb                *mygrid;
   Electron_Operators  *myelectron;

   double *U, *Vt, *S;

public:

   /* Constructors */
   Geodesic(int minimizer0, Molecule& mymolecule) {
      minimizer  = minimizer0;
      myelectron = mymolecule.myelectron;
      mygrid     = mymolecule.mygrid;
      U  = mygrid->g_allocate(1);
      Vt = mygrid->m_allocate(-1,1);
      S = new double[mygrid->ne[0]+mygrid->ne[1]];
   }


   /* destructor */
   ~Geodesic() {
         delete [] U;
         delete [] Vt;
         delete [] S;
    }


   double start(double *A, double *max_sigma) {
      double *V = mygrid->m_allocate(-1,1);
      mygrid->ggm_SVD(A,U,S,V);

      int neall = mygrid->ne[0] + mygrid->ne[1];
      double msig =0.0;
      for (int i=0; i<neall; ++i)
         if (fabs(S[i])>msig)
            msig = fabs(S[i]);
      *max_sigma = msig;

      /* calculate Vt */
      mygrid->mm_transpose(-1,V,Vt);

      delete [] V;

      /* calculate  and return 2*<A|H|psi> */
      return(2.0*myelectron->eorbit(A));

    }
 

};

#endif

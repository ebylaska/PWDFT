#ifndef	_COULOMB_H_
#define _COULOMB_H_




#include	"Pneb.hpp"

class	Coulomb_Operator {

   double *vg;
   Pneb   *mypneb;

public:

   /* Constructors */
   Coulomb_Operator(Pneb *);

   /* destructor */
   ~Coulomb_Operator() {
         delete [] vg;
    }

    void   vcoulomb(double *, double *);
    double ecoulomb(double *);
};

#endif

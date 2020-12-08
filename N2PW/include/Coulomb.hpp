#ifndef	_COULOMB_H_
#define _COULOMB_H_


using namespace std;

#include	"Pneb.hpp"

class	Coulomb_Operator {

   float *vg;
   Pneb   *mypneb;

public:

   /* Constructors */
   Coulomb_Operator(Pneb *);

   /* destructor */
   ~Coulomb_Operator() {
         delete [] vg;
    }

    void   vcoulomb(float *, float *);
    float ecoulomb(float *);
};

#endif

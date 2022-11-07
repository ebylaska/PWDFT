#ifndef	_KINETIC_H_
#define _KINETIC_H_




#include	"Pneb.hpp"

class	Kinetic_Operator {

   double *tg;
   Pneb   *mypneb;

public:

   /* Constructors */
   Kinetic_Operator(Pneb *);

   /* destructor */
   ~Kinetic_Operator() {
         delete [] tg;
    }

    void   ke(double *, double *);
    double ke_ave(double *);
};

#endif

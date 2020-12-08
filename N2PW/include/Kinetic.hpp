#ifndef	_KINETIC_H_
#define _KINETIC_H_


using namespace std;

#include	"Pneb.hpp"

class	Kinetic_Operator {

   float *tg;
   Pneb   *mypneb;

public:

   /* Constructors */
   Kinetic_Operator(Pneb *);

   /* destructor */
   ~Kinetic_Operator() {
         delete [] tg;
    }

    void   ke(float *, float *);
    float ke_ave(float *);
};

#endif

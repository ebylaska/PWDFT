#ifndef	_IONS_H_
#define _IONS_H_


#include        <string>
using namespace std;

#include	"rtdb.hpp"

class	Ion {

   char   *atomarray;

public:
   int nion,nkatm;
   int *katm,*natm;
   double *charge;
   double *mass;
   double *rion1,*rion2;

   /* Constructors */
   Ion(RTDB&);
   Ion(string);

   /* destructor */
   ~Ion() {
      delete [] katm;
      delete [] natm;
      delete [] atomarray;
      delete [] charge;
      delete [] mass;
      delete [] rion1;
      delete [] rion2;
    }

    char *symbol(const int i) { return &atomarray[3*katm[i]]; }
    char *atom(const int ia)  { return &atomarray[3*ia]; }

};

#endif

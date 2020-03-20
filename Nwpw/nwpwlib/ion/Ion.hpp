#ifndef	_IONS_H_
#define _IONS_H_


#include        <string>

#include	"rtdb.hpp"

class	Ion {

   char   *atomarray;

public:
   int nion,nkatm;
   int *katm,*natm;
   double *charge;
   double *mass;
   double *rion0,*rion1,*rion2;

   /* Constructors */
   Ion(RTDB&);
   Ion(std::string);

   /* destructor */
   ~Ion() {
      delete [] katm;
      delete [] natm;
      delete [] atomarray;
      delete [] charge;
      delete [] mass;
      delete [] rion0;
      delete [] rion1;
      delete [] rion2;
    }

    void shift() 
    { 
       for (auto i=0; i<(3*nion); ++i) rion0[i] = rion1[i];
       for (auto i=0; i<(3*nion); ++i) rion1[i] = rion2[i];
    }
    char *symbol(const int i) { return &atomarray[3*katm[i]]; }
    char *atom(const int ia)  { return &atomarray[3*ia]; }

};

#endif

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
   double *dti;
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
      delete [] dti;
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

    void optimize_step(const double *fion) 
    {
       for (auto ii=0; ii<nion; ++ii)
       {
          double alpha = sqrt(dti[ii]);
          rion2[3*ii]   = rion1[3*ii]   + alpha*fion[3*ii];
          rion2[3*ii+1] = rion1[3*ii+1] + alpha*fion[3*ii+1];
          rion2[3*ii+2] = rion1[3*ii+2] + alpha*fion[3*ii+2];
       }
    }

};

#endif

#ifndef	_IONS_H_
#define _IONS_H_


#include        <string>
#include        <cmath>

#include	"rtdb.hpp"
#include	"Control2.hpp"

class	Ion {

   char   *atomarray;

public:
   int ion_total_charge;
   int nion,nkatm; // number of ions
   int *katm; // element ID number
   int *natm;
   double *charge,*zv_psp;
   double *mass;
   double *dti;
   double *rion0,*rion1,*rion2; // coordinates of ions

   /* Constructors */
   Ion(RTDB&, Control2&);
   Ion(std::string, Control2&);

   /* destructor */
   ~Ion() {
      delete [] katm;
      delete [] natm;
      delete [] atomarray;
      delete [] charge;
      delete [] zv_psp;
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
    double amu(const int i) { return mass[i]/1822.89; }
    void writejsonstr(string&);
    double rion(int i, int ii) { return rion1[3*ii+i];}
    double vion(int i, int ii) { return rion0[3*ii+i];}
    double gc(int i) 
    { 
       double ss = 0.0;
       for (auto ii=0; ii<nion; ++ii) 
          ss += rion1[3*ii+i];
       return ss/((double) nion);
    }
    double vgc(int i) 
    { 
       double ss = 0.0;
       for (auto ii=0; ii<nion; ++ii) 
          ss += rion0[3*ii+i];
       return ss/((double) nion);
    }
    double com(int i) 
    {
       double tmass = 0.0;
       double ss    = 0.0;
       for (auto ii=0; ii<nion; ++ii)
       { 
          tmass += mass[ii];
          ss    += mass[ii]*rion1[3*ii+i];
       }
       return ss/tmass;
    }
    double vcom(int i) 
    {
       double tmass = 0.0;
       double ss    = 0.0;
       for (auto ii=0; ii<nion; ++ii)
       { 
          tmass += mass[ii];
          ss    += mass[ii]*rion0[3*ii+i];
       }
       return ss/tmass;
    }

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

    void set_zv_psp(const int ia, const double zv) { zv_psp[ia] = zv; }
    double total_zv()
    {
       double tcharge = 0.0;
       for (auto ia=0; ia< nkatm; ++ia)
          tcharge += zv_psp[ia]*natm[ia];
       return tcharge;
    }

};

#endif

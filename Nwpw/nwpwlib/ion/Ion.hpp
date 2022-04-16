#ifndef	_IONS_H_
#define _IONS_H_


#include        <string>
#include        <cmath>

#include	"rtdb.hpp"
#include	"Control2.hpp"
#include	"ion_rcovalent.hpp"

namespace pwdft {
using namespace pwdft;

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


    void fixed_step(const double alpha, const double *fion) 
    {
       for (auto i=0; i<(3*nion); ++i)
          rion2[i] = rion1[i] + alpha*fion[i];
    }


    void set_zv_psp(const int ia, const double zv) { zv_psp[ia] = zv; }
    double total_zv()
    {
       double tcharge = 0.0;
       for (auto ia=0; ia< nkatm; ++ia)
          tcharge += zv_psp[ia]*natm[ia];
       return tcharge;
    }

    double xrms()
    {
       double dx,dy,dz;
       double xx=0.0;
       for (auto ii=0; ii<nion; ++ii)
       {
          dx = rion1[3*ii]   - rion0[3*ii];
          dy = rion1[3*ii+1] - rion0[3*ii+1];
          dz = rion1[3*ii+2] - rion0[3*ii+2];
          xx += dx*dx + dy*dy + dz*dz;
       }
       return sqrt(xx)/((double) nion);
    }

    double xmax()
    {
       double dx,dy,dz,x,xx;
       double yy = 0.0;
       for (auto ii=0; ii<nion; ++ii)
       {
          dx = rion1[3*ii]   - rion0[3*ii];
          dy = rion1[3*ii+1] - rion0[3*ii+1];
          dz = rion1[3*ii+2] - rion0[3*ii+2];
          xx = dx*dx + dy*dy + dz*dz;
          x = sqrt(xx); if (x>yy) yy = x;
       }
       return yy;
    }

    std::string print_bond_angle_torsions() { return ion_print_bond_angle_torsions(nion,katm,atomarray,rion1); }


    double com_fion(const double *fion, const int i) {
       double tmass = 0.0;
       double ss    = 0.0;
       for (auto ii=0; ii<nion; ++ii) { 
          tmass += mass[ii];
          ss    += mass[ii]*fion[3*ii+i];
       }
       return ss/tmass;
    }

    void remove_com_fion(double *fion) {
       double tmass = 0.0;
       double gx = 0.0;
       double gy = 0.0;
       double gz = 0.0;
       for (auto ii=0; ii<nion; ++ii) { 
          tmass += mass[ii];
          gx    += mass[ii]*fion[3*ii];
          gy    += mass[ii]*fion[3*ii+1];
          gz    += mass[ii]*fion[3*ii+2];
       }
       gx /= tmass;
       gy /= tmass;
       gz /= tmass;
       for (auto ii=0; ii<nion; ++ii) { 
          fion[3*ii]   -= gx;
          fion[3*ii+1] -= gy;
          fion[3*ii+2] -= gz;
       }
    }

    double max_fion(const double *fion) {
       double x,xx;
       double y = 0.0;
       for (auto ii=0; ii<nion; ++ii) { 
          xx = fion[3*ii]*fion[3*ii] + fion[3*ii+1]*fion[3*ii+1] + fion[3*ii+2]*fion[3*ii+2];
          x = sqrt(xx); if (x>y) y = x;
       }
       return y;
   }


    double rms_fion(const double *fion) {
       double ss    = 0.0;
       for (auto ii=0; ii<nion; ++ii) { 
          ss += fion[3*ii]*fion[3*ii] + fion[3*ii+1]*fion[3*ii+1] + fion[3*ii+2]*fion[3*ii+2];
       }
       return sqrt(ss)/((double) nion);
   }



};
}

#endif

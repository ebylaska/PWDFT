#ifndef	_ELECTRON_HPP_
#define _ELECTRON_HPP_


using namespace std;

#include	"Pneb.hpp"
#include	"Kinetic.hpp"
#include	"Coulomb.hpp"
#include        "Strfac.hpp"
#include	"Pseudopotential.hpp"

class	Electron_Operators {

   Pneb   *mygrid;

   Kinetic_Operator *myke;
   Coulomb_Operator *mycoulomb;
   Pseudopotential  *mypsp;

   double *Hpsi, *psi_r, *vl, *vall, *vc, *xcp, *xce, *x, *hmltmp;

   double omega,scal2,scal1,dv;

   int ispin,neall,n2ft3d,shift1,shift2;
   int counter=0;

public:

   /* Constructors */
   Electron_Operators(Pneb *, Kinetic_Operator *, Coulomb_Operator *, Pseudopotential *);

   /* destructor */
   ~Electron_Operators() {
         delete [] Hpsi;
         delete [] psi_r;
         delete [] vl;
         delete [] x;
         delete [] vall;
         delete [] vc;
         delete [] xcp;
         delete [] xce;
         delete [] hmltmp;
    }

    void genrho(double *, double *);
    void run(double *, double *, double *, double *);
    void get_Tgradient(double *, double *, double *);
    void get_Gradient(double *);
    void gen_Tangent(double *, double *, double *);
    void gen_hml(double *, double *);
    void gen_Hpsi_k(double *);
    void gen_psi_r(double *);
    void gen_density(double *);
    void gen_densities(double *, double *, double *);
    void gen_scf_potentials(double *, double *, double *);
    void gen_vl_potential();

    double psi_vl_ave(double *);
    double psi_vnl_ave(double *);

    double energy(double *, double *, double *, double *);
    double eorbit(double *);
    double ehartree(double *);
    double exc(double *);
    double pxc(double *);

    //void semicoreforce(double *);
    
};

#endif
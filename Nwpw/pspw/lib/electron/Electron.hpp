#pragma once

#include	"Pneb.hpp"
#include	"Kinetic.hpp"
#include	"Coulomb.hpp"
#include	"exchange_correlation.hpp"
#include        "Strfac.hpp"
#include	"Pseudopotential.hpp"

namespace pwdft {


class	Electron_Operators {

   Pneb   *mygrid;

   Kinetic_Operator *myke;
   Coulomb_Operator *mycoulomb;
   XC_Operator      *myxc;
   Pseudopotential  *mypsp;

   double *Hpsi, *psi_r, *vl, *vall, *vc, *xcp, *xce, *x, *hmltmp;

   double omega,scal2,scal1,dv;

   int ispin,neall,n2ft3d,shift1,shift2,npack1;

public:
   int counter=0;

   /* Constructors */
   Electron_Operators(Pneb *, Kinetic_Operator *, Coulomb_Operator *, XC_Operator *, Pseudopotential *);

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
    void semicore_density_update();

    double vl_ave(double *);
    double vnl_ave(double *);

    double energy(double *, double *, double *, double *);
    double eorbit(double *);
    double ehartree(double *);
    double exc(double *);
    double pxc(double *);
    double eke(double *);

    void gen_energies_en(double *, double *, double *, double *, double *, double *);

    void add_dteHpsi(double, double *, double *);

    void vl_force(double *, double *);
    void vnl_force(double *, double *);
    void semicore_force(double *);

    void apc_force(double *, double *);

    bool is_v_apc_on() { return mypsp->myapc->v_apc_on; }
    
};

}

#ifndef _ELECTRON_HPP_
#define _ELECTRON_HPP_

#pragma once

#include	"Pneb.hpp"
#include	"Kinetic.hpp"
#include	"Coulomb12.hpp"
#include	"exchange_correlation.hpp"
#include        "Strfac.hpp"
#include	"Pseudopotential.hpp"

namespace pwdft {


class	Electron_Operators {

   Pneb   *mygrid;

   Kinetic_Operator   *myke;
   Coulomb12_Operator *mycoulomb12;
   XC_Operator        *myxc;
   Pseudopotential    *mypsp;

   double *Hpsi, *psi_r, *vl, *vall, *vc, *xcp, *xce, *x, *rho, *hmltmp;
   double *vlr_l;

   double omega,scal2,scal1,dv;

   int ispin,neall,n2ft3d,shift1,shift2,npack1;
   bool aperiodic = false;
   bool  periodic = false;

public:
   int counter=0;

   /* Constructors */
   Electron_Operators(Pneb *, Kinetic_Operator *, Coulomb12_Operator *, XC_Operator *, Pseudopotential *);

   /* destructor */
   ~Electron_Operators() {
         delete [] Hpsi;
         delete [] psi_r;
         delete [] vl;
         delete [] x;
         delete [] rho;
         delete [] vall;
         delete [] vc;
         delete [] xcp;
         delete [] xce;
         delete [] hmltmp;
         if (aperiodic) delete [] vlr_l;
    }

    void genrho(double *, double *);
    void run(double *, double *, double *, double *);
    void get_Tgradient(double *, double *, double *);
    void get_Gradient(double *);
    void gen_Tangent(double *, double *, double *);
    void gen_hml(double *, double *);
    void gen_hmlt(double *, double *);
    void gen_Hpsi_k(double *);
    void gen_psi_r(double *);
    void gen_density(double *);
    void gen_densities(double *, double *, double *);
    void gen_scf_potentials(double *, double *, double *);
    void gen_vl_potential();
    void semicore_density_update();

    double vl_ave(double *);
    double vlr_ave(double *);
    double vnl_ave(double *);

    double energy(double *, double *, double *, double *);
    double eorbit(double *);
    double ehartree(double *);
    double ehartree2(double *);
    double exc(double *);
    double pxc(double *);
    double eke(double *);

    void gen_energies_en(double *, double *, double *, double *, double *, double *);

    void add_dteHpsi(double, double *, double *);

    void vl_force(double *, double *);
    void vlr_force(double *, double *);
    void vnl_force(double *, double *);
    void semicore_force(double *);

    void apc_force(double *, double *);

    bool is_v_apc_on() { return mypsp->myapc->v_apc_on; }
    bool is_aperiodic() { return aperiodic; }
    bool is_periodic() { return periodic; }
    
};

}

#endif

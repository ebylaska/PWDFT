
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Control2.hpp"
#include "Geodesic.hpp"
#include "Ion.hpp"
#include "Molecule.hpp"
#include "Parallel.hpp"
#include "Pneb.hpp"
#include "util_date.hpp"
#include "util_linesearch.hpp"
#include "nwpw_scf_mixing.hpp"

namespace pwdft {

/* create dummy function call to Geodesic class functions */
static Geodesic *mygeodesic_ptr;
static double dummy_energy(double t) { return mygeodesic_ptr->energy(t); }
static double dummy_denergy(double t) { return mygeodesic_ptr->denergy(t); }

/******************************************
 *                                        *
 *            cgsd_bybminimize2           *
 *                                        *
 ******************************************/
double cgsd_bybminimize2(Molecule &mymolecule, Geodesic *mygeodesic, double *E,
                        double *deltae, double *deltac, int current_iteration,
                        int ks_it_in, int ks_it_out, 
                        nwpw_scf_mixing &scfmix,
                        double tole, double tolc) 
{
   bool done = false;
   double tmin = 0.0;
   double deltat_min = 1.0e-3;
   double deltat;
   double sum0, sum1, scale, total_energy;
   double dE, max_sigma, min_sigma;
   double Eold, dEold, Enew;
   double tmin0, deltae0, perror,error_out;
 
   Pneb *mygrid = mymolecule.mygrid;
   mygeodesic_ptr = mygeodesic;

   double ks_deltae = tole;
   int ispin = mygrid->ispin;
   int *neq = mygrid->neq;
   double *vall = mygrid->r_nalloc(ispin);
   mymolecule.gen_vall();
   mymolecule.get_vall(vall);



   // ion-ion energy 
   //double eion = mymolecule.eion();

   //**********************
   //**** bybminimizer ****
   //**********************
   double e0 = mymolecule.psi_KS_update(ks_it_in,ks_deltae,perror,vall,
                                        ispin, neq, mymolecule.psi1,
                                        deltae, std::cout);

   /* iniitialize blocked cg */
 
   // Making an extra call to electron.run and energy
   total_energy = mymolecule.gen_all_energies();
 
   //|-\____|\/-----\/\/->    End Parallel Section    <-\/\/-----\/|____/-|
 
   //mygrid->g_deallocate(H0);
   //mygrid->g_deallocate(G1);
 
   return total_energy;
}

} // namespace pwdft

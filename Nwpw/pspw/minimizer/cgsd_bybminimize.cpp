
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
 *            cgsd_bybminimize            *
 *                                        *
 ******************************************/
double cgsd_bybminimize(Molecule &mymolecule, Geodesic *mygeodesic, double *E,
                        double *deltae, double *deltac, int current_iteration,
                        int it_in, double tole, double tolc) 
{
   bool done = false;
   double tmin = 0.0;
   double deltat_min = 1.0e-3;
   double deltat;
   double sum0, sum1, scale, total_energy;
   double dE, max_sigma, min_sigma;
   double Eold, dEold, Enew;
   double tmin0, deltae0;
 
   Pneb *mygrid = mymolecule.mygrid;
   mygeodesic_ptr = mygeodesic;
   double *vall_out;

   //int sd_it = 10
   //int cg_it = 10
 //   if (set_iterations) then
 //       it_in = iterations
 //       sd_it = 2 
 //       cg_it = 1 
 //     else
 //       it_in = control_it_in()*control_it_out()
 //       sd_it = 10
 //       cg_it = 10
 //     end if

//      int maxit_orbs = control_ks_maxit_orbs() 
      //precondition = control_precondition()
      //ispin = control_ispin()
      double deltav_old = 10.0;
      double deltav     = 0.0;
 
 //     stalled       = false;
  //    std::vector<double> deltae_history;
   //   deltae_history.push_back(0.0);
  //    deltae_history.push_back(0.0);
   //   deltae_history.push_back(0.0);
    //  deltae_history.push_back(0.0);
      int stalled_count     = 0;
      int sd_count          = 0;

 //     vall_in = 
  //    vall_out = 
 //     vall_junk = 
 //     rho_in = 

   //**********************
   //**** bybminimizer ****
   //**********************

   nwpw_scf_mixing scfmix(mygrid->d3db::parall,1,0.10,2,mygrid->ispin,mygrid->n2ft3d,vall_out);


 
   /* iniitialize blocked cg */
 
   // Making an extra call to electron.run and energy
   total_energy = mymolecule.gen_all_energies();
 
   //|-\____|\/-----\/\/->    End Parallel Section    <-\/\/-----\/|____/-|
 
   //mygrid->g_deallocate(H0);
   //mygrid->g_deallocate(G1);
 
   return total_energy;
}

} // namespace pwdft

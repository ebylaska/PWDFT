
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Control2.hpp"
#include "band_Geodesic.hpp"
#include "Ion.hpp"
#include "Solid.hpp"
#include "Parallel.hpp"
#include "Cneb.hpp"
#include "util_date.hpp"
#include "util_linesearch.hpp"

namespace pwdft {

/* create dummy function call to Geodesic class functions */
static band_Geodesic *mygeodesic_ptr;
static double dummy_energy(double t) { return mygeodesic_ptr->energy(t); }
static double dummy_denergy(double t) { return mygeodesic_ptr->denergy(t); }

/******************************************
 *                                        *
 *            band_cgsd_bybminimize2      *
 *                                        *
 ******************************************/
double band_cgsd_bybminimize2(Solid &mysolid, band_Geodesic *mygeodesic, double *E,
                              double *deltae, double *deltac, int current_iteration,
                              int ks_it_in, int ks_it_out,
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
 
   Cneb *mygrid = mysolid.mygrid;
   mygeodesic_ptr = mygeodesic;

   double ks_deltae = tole;
   int ispin = mygrid->ispin;
   int *neq = mygrid->neq;
   int nbrillq = mygrid->nbrillq;
   double *vall = mygrid->c_nalloc(ispin);
   mysolid.gen_vall();
   mysolid.get_vall(vall);


   //**********************
   //**** bybminimizer ****
   //**********************
   double e0 = mysolid.cpsi_KS_update(ks_it_in,ks_deltae,perror,vall,
                                      ispin, neq, nbrillq, mysolid.psi1,
                                      deltae, std::cout);

   /* iniitialize blocked cg */
 
   // Making an extra call to electron.run and energy
   total_energy = mysolid.gen_all_energies();
 
   //|-\____|\/-----\/\/->    End Parallel Section    <-\/\/-----\/|____/-|
 
 
   return total_energy;
}


} // namespace pwdft

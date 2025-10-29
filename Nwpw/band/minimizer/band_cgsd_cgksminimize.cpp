
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
static double dummy_energy(double t) { return mygeodesic_ptr->energy0(t); }    // needs to be modified
static double dummy_denergy(double t) { return mygeodesic_ptr->denergy(t); }  // needs to be modified

/******************************************
 *                                        *
 *            band_cgsd_cgksminimize      *
 *                                        *
 ******************************************/
/* Only minimize the orbitals, keeping the density and ks potetnials fixed
*/
double band_cgsd_cgksminimize(Solid &mysolid, band_Geodesic *mygeodesic, double *E,
                              double *deltae, double *deltac, int current_iteration,
                              int it_in, double tole, double tolc) {
  bool done = false;
  double tmin = 0.0;
  double deltat_min = 1.0e-3;
  double deltat;
  double sum0, sum1, scale, total_energy;
  double dE, max_sigma, min_sigma;
  double Eold, dEold, Enew;
  double tmin0, deltae0;

  Cneb *mygrid = mysolid.mygrid;
  mygeodesic_ptr = mygeodesic;


  /* get the initial gradient and direction */
  double *G1 = mygrid->g_allocate_nbrillq_all();
  double *H0 = mygrid->g_allocate_nbrillq_all();

  //|-\____|\/-----\/\/->    Start Parallel Section    <-\/\/-----\/|____/-|

  total_energy = mysolid.psi_1get_Tgradient0(G1); // needs to be modified
  sum1 = mygrid->gg_traceall(G1, G1);
  Enew = total_energy;
  *deltac = sum1;

  mygrid->gg_copy(G1, H0);

  /******************************************
   ****                                  ****
   **** Start of conjugate gradient loop ****
   ****                                  ****
   ******************************************/
  int it = 0;
  tmin = deltat_min;
  while ((!done) && ((it++) < it_in)) {
    /* initialize the geoedesic line data structure */
    dEold = mygeodesic->start(H0, &max_sigma, &min_sigma);


    /* line search */
    if (tmin > deltat_min)
      deltat = tmin;
    else
      deltat = deltat_min;

    tmin0 = tmin;
    deltae0 = *deltae;

    Eold = Enew;

    Enew = util_linesearch(0.0, Eold, dEold, deltat, &dummy_energy,
                           &dummy_denergy, 0.50, &tmin0, &deltae0, 2);
    tmin = tmin0;
    *deltae = deltae0;
    
    mygeodesic->psi_final(tmin);

    /* exit loop early */
    done = ((it >= it_in) || ((std::fabs(*deltae) < tole) && (*deltac < tolc)));

    /* transport the previous search directions */
    mygeodesic->psi_1transport(tmin, H0);

    /* make psi1 <--- psi2(tmin) */
    mysolid.swap_psi1_psi2();

    if (!done) {
      /* get the new gradient - but do not updates densities */
      total_energy = mysolid.psi_1get_Tgradient0(G1); // needs to be modified
      sum0 = sum1;
      sum1 = mygrid->gg_traceall(G1, G1);
      *deltac = sum1;

      /* the new direction using Fletcher-Reeves */
      if ((std::fabs(*deltae) <= (1.0e-2)) && (tmin > deltat_min)) {
        if (sum0 > 1.0e-9)
          scale = sum1 / sum0;
        else
          scale = 0.0;

        mygrid->g_Scale(scale, H0);
        mygrid->gg_Sum2(G1, H0);
      }

      /* the new direction using steepest-descent */
      else
        mygrid->gg_copy(G1, H0);

      // mygrid->gg_copy(G1,H0);
    }
  }
  // Making an extra call to electron.run and energy
  total_energy = mysolid.energy0();

  //|-\____|\/-----\/\/->    End Parallel Section    <-\/\/-----\/|____/-|

  mygrid->g_deallocate(H0);
  mygrid->g_deallocate(G1);

  return total_energy;
}

} // namespace pwdft

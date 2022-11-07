
#include        <iostream>
#include        <cstdio>
#include        <cmath>
#include        <cstdlib>
#include        <string>

#include        "util_linesearch.hpp"
#include        "Parallel.hpp"
#include        "Control2.hpp"
#include        "util_date.hpp"
#include        "Ion.hpp"
#include        "Pneb.hpp"
#include        "Molecule.hpp"
#include        "Geodesic.hpp"
#include        "pspw_lmbfgs.hpp"


namespace pwdft {

/* create dummy function call to Geodesic class functions */
static Geodesic *mygeodesic_ptr;
static double dummy_energy(double t)  { return mygeodesic_ptr->energy(t); }
static double dummy_denergy(double t) { return mygeodesic_ptr->denergy(t); }


/******************************************
 *                                        *
 *            cgsd_bfgsminimize           *
 *                                        *
 ******************************************/
double cgsd_bfgsminimize(Molecule& mymolecule, Geodesic& mygeodesic, pspw_lmbfgs &psi_lmbfgs, 
                         double *E, double *deltae, double *deltac, 
                         int current_iteration, int it_in, double tole, double tolc)
{
   bool   done = false;
   double tmin = 0.0;
   double deltat_min=1.0e-2;
   double deltat;
   double sum0,sum1,scale,total_energy;
   double dE,max_sigma,min_sigma;
   double Eold,dEold,Enew;
   double tmin0,deltae0;

   Pneb *mygrid = mymolecule.mygrid;
   mygeodesic_ptr = &mygeodesic;



   /* get the initial gradient and direction */
   double *G0 = mygrid->g_allocate(1);
   double *S0 = mygrid->g_allocate(1);

   total_energy  = mymolecule.psi_1get_Tgradient(G0);
   sum1 = mygrid->gg_traceall(G0,G0);
   Enew = total_energy;

   if (current_iteration==0) 
   {
      psi_lmbfgs.start(G0);
      mygrid->gg_copy(G0,S0);
   }
   else
   {
      psi_lmbfgs.fetch(tmin,G0,S0);

      // reset to gradient if <S0|G0> <= 0.0
      double kappa = mygrid->gg_traceall(G0,S0);
      if (kappa<=0.0) mygrid->gg_copy(G0,S0);
   }


   /******************************************
    ****                                  ****
    ****   Start of BFGS iteration loop   ****
    ****                                  ****
    ******************************************/
   int it = 0;
   tmin  = deltat_min;
   while ((!done) && ((it++) < it_in))
   {
      /* initialize the geoedesic line data structure */
      dEold = mygeodesic.start(S0,&max_sigma, &min_sigma);

      /* line search */
      if (tmin > deltat_min)
         deltat = tmin;
      else
         deltat = deltat_min;

      tmin0   = tmin;
      deltae0 = *deltae;

      Eold = Enew;
      Enew = util_linesearch(0.0,Eold,dEold,deltat,
                       &dummy_energy, &dummy_denergy,
                       0.50,&tmin0,&deltae0,2);
      tmin = tmin0;
      *deltae = deltae0;
      *deltac = mymolecule.rho_error();
      mygeodesic.psi_final(tmin);


      /* exit loop early */
      done = ((it >= it_in) || ((std::fabs(*deltae)<tole) && (*deltac<tolc)));

      /* make psi1 <--- psi2(tmin) */
      mymolecule.swap_psi1_psi2();

      if (!done)
      {
         /* get the new gradient - also updates densities */
         total_energy  = mymolecule.psi_1get_Tgradient(G0);
         psi_lmbfgs.fetch(tmin,G0,S0);

         // reset to gradient if <S0|G0> <= 0.0
         double kappa = mygrid->gg_traceall(G0,S0);
         if (kappa<=0.0) mygrid->gg_copy(G0,S0);

      }

   }
   // Making an extra call to electron.run and energy
   total_energy  = mymolecule.gen_all_energies();

   mygrid->g_deallocate(S0);
   mygrid->g_deallocate(G0);
   
   return total_energy;
}

}

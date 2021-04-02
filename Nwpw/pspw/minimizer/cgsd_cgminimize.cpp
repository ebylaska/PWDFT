
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


static Geodesic *mygeodesic_ptr;
static double dummy_energy(double t) { return mygeodesic_ptr->energy(t); }
static double dummy_denergy(double t) { return mygeodesic_ptr->denergy(t); }


double cgsd_cgminimize(Molecule& mymolecule, Geodesic& mygeodesic, double *E, double *deltae, double *deltac, 
                       int current_iteration, int it_in, double tole, double tolc)
{
   bool   done = false;
   double tmin = 0.0;
   double deltat_min=1.0e-3;
   double deltat;
   double sum0,sum1,scale,total_energy;
   double dE,max_sigma,min_sigma;
   double Eold,dEold,Enew;
   double tmin0,deltae0;

   Pneb *mygrid = mymolecule.mygrid;
   mygeodesic_ptr = &mygeodesic;


   /* get the initial gradient and direction */
   double *G1 = mygrid->g_allocate(1);
   double *H0 = mygrid->g_allocate(1);

   total_energy  = mymolecule.psi_1get_Tgradient(G1);
   sum1 = mygrid->gg_traceall(G1,G1);
   Enew = total_energy;

   //std::cout << "SUM1=" << sum1 << " ENERGY = " << total_energy << std::endl;

   mygrid->gg_copy(G1,H0);


   /******************************************
    ****                                  ****
    **** Start of conjugate gradient loop ****
    ****                                  ****
    ******************************************/
   int it = 0;
   tmin  = deltat_min;
   while ((!done) && ((it++) < it_in))
   {
      /* initialize the geoedesic line data structure */
      dEold = mygeodesic.start(H0,&max_sigma, &min_sigma);
      //std::cout << "MAX_SIGMA=" << max_sigma << " MIN_SIGMA="<< min_sigma << std::endl;

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

      //std::cout << "out geodesic.energy Eold=" << Eold << " dEold=" << dEold << std::endl;
      //std::cout << "out geodesic.energy Enew=" << Enew << " DELTAE = " << *deltae << " tmin=" << tmin << std::endl;

      /* exit loop early */
      done = ((it >= it_in) || ((fabs(*deltae)<tole) && (*deltac<tolc)));


      //std::cout << "it=" << it << " done=" << done << std::endl;
      if (!done)
      {
         /* transport the previous search directions */
         mymolecule.swap_psi1_psi2();

         /* make psi1 <--- psi2(tmin) */

         /* get the new gradient - also updates densities */
         total_energy  = mymolecule.psi_1get_Tgradient(G1);
         sum0 = sum1;
         sum1 = mygrid->gg_traceall(G1,G1);
         //std::cout << "SUM0=" << sum0 <<  " SUM1=" << sum1 << std::endl;

         /* the new direction using Fletcher-Reeves */
         if ( (fabs(*deltae)<=(1.0e-2)) &&  (tmin>deltat_min))
         {
             if (sum0>1.0e-9)
                scale = sum1/sum0;
             else
                scale = 0.0;

             mygrid->g_Scale(scale,H0);
             mygrid->gg_Sum2(G1,H0);
         } 

         /* the new direction using steepest-descent */
         else
              mygrid->gg_copy(G1,H0);

         //mygrid->gg_copy(G1,H0);
      }

   }

   total_energy  = mymolecule.gen_all_energies();


   mygrid->g_deallocate(H0);
   mygrid->g_deallocate(G1);
   
   return total_energy;
}

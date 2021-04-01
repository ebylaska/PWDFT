
#include        <iostream>
#include        <cstdio>
#include        <cmath>
#include        <cstdlib>
#include        <string>

#include        "Parallel.hpp"
#include        "Control2.hpp"
#include        "util_date.hpp"
#include        "Ion.hpp"
#include        "Pneb.hpp"
#include        "Molecule.hpp"
#include        "Geodesic.hpp"

double cgsd_cgminimize(Molecule& mymolecule, Geodesic& mygeodesic, double *E, double *deltae, double *deltac, int current_iteration, int it_in)
{
   bool   done = false;
   double tmin = 0.0;;
   double deltat_min=1.0e-3;
   double sum0,sum1,scale,total_energy;
   double dE,max_sigma;

   Pneb *mygrid = mymolecule.mygrid;


   /* get the initial gradient and direction */
   double *G1 = mygrid->g_allocate(1);
   double *H0 = mygrid->g_allocate(1);

   total_energy  = mymolecule.psi_1get_Tgradient(G1);
   sum1 = mygrid->gg_traceall(G1,G1);

   std::cout << "SUM1=" << sum1 << std::endl;

   mygrid->gg_copy(G1,H0);


   /******************************************
    ****                                  ****
    **** Start of conjugate gradient loop ****
    ****                                  ****
    ******************************************/
   int it = 0;
   while ((!done) && ((it++) < it_in))
   {
      /* initialize the geoedesic line data structure */
      //dE = mygeodesic.start(H0,&max_sigma);

      /* line search */

      /* exit loop early */
      done = (it >= it_in);

      std::cout << "it=" << it << " done=" << done << std::endl;
      if (!done)
      {
         /* transport the previous search directions */

         /* make psi1 <--- psi2(tmin) */

         /* get the new gradient - also updates densities */
         total_energy  = mymolecule.psi_1get_Tgradient(G1);
         sum0 = sum1;
         sum1 = mygrid->gg_traceall(G1,G1);
         std::cout << "SUM0=" << sum0 <<  " SUM1=" << sum1 << std::endl;

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
      }

   }

   total_energy  = mymolecule.gen_all_energies();

   *deltae = 1.33434340e-4;
   *deltac = 2.33434211e-6;

   mygrid->g_deallocate(H0);
   mygrid->g_deallocate(G1);
   
   return total_energy;
}

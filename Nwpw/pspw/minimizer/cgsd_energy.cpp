#include        <iostream>
#include        <cstdio>
#include        <cmath>
#include        <cstdlib>
#include        <string>

#include        "Parallel.hpp"
#include        "Control2.hpp"
#include        "Ion.hpp"
#include        "Pneb.hpp"
#include        "Molecule.hpp"

/******************************************
 *                                        *
 *            cgsd_noit_energy            *
 *                                        *
 ******************************************/
double cgsd_noit_energy(Molecule& mymolecule) 
{
   Parallel *parall = mymolecule.mygrid->d3db::parall;

   double total_energy  = mymolecule.gen_all_energies();

   /* report summary of results */
   if (parall->is_master())
   {
      std::cout << "     ==================  optimization turned off  ===================\n" << std::endl;
      std::cout << std::endl;
      std::cout << mymolecule;
   }

   return total_energy;
}


/******************************************
 *                                        *
 *              cgsd_energy               *
 *                                        *
 ******************************************/
double cgsd_energy(Control2& control, Molecule& mymolecule)
{
   Parallel *parall = mymolecule.mygrid->d3db::parall;
   Pneb *mygrid = mymolecule.mygrid;
   Ion  *myion  = mymolecule.myion;

   int    ne[2],ispin,nion;
   double E[50],total_energy;

   int it_in   = control.loop(0);
   int it_out  = control.loop(1);
   double tole = control.tolerances(0);
   double tolc = control.tolerances(1);

   int minimizer = control.minimizer();

   /* generate phaze factors and local psp and core density */
/*
      call phafac()
      call  ewald_phafac()
      call electron_gen_vl_potential()
      if (psp_semicore(0)) call semicore_density_update()
*/

   std::cout << "cgsd_energy: minimizer = " << minimizer << std::endl;

   total_energy = 0.0;

   return total_energy;
}


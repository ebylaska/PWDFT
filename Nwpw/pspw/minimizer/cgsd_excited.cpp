#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Control2.hpp"
#include "Geodesic12.hpp"
#include "Ion.hpp"
#include "Molecule.hpp"
#include "Parallel.hpp"
#include "Cneb.hpp"
#include "iofmt.hpp"
//#include "band_lmbfgs.hpp"
//#include "band_lmbfgs2.hpp"
#include "util_date.hpp"

#include "cgsd_excited.hpp"

#include "iofmt.hpp"

namespace pwdft {


/******************************************
 *                                        *
 *              cgsd_excited              *
 *                                        *
 ******************************************/
void cgsd_excited(Control2 &control, Molecule &mymolecule, bool doprint, std::ostream &coutput) 
{
   Parallel *parall = mymolecule.mygrid->d3db::parall;
   Pneb *mygrid = mymolecule.mygrid;
   Ion *myion = mymolecule.myion;
 
   bool stalled = false;
   double E[70],total_energy,deltae,deltae_old,deltac;
 
   int it_in   = control.loop(0);
   int it_out  = control.loop(1);
   double tole = control.tolerances(0);
   double tolc = control.tolerances(1);
 
   double dt = control.time_step();
   double dte = dt/sqrt(control.fake_mass());
 
   int minimizer   = control.minimizer();
   int lmbfgs_size = control.lmbfgs_size();
 
   bool hprint = (parall->is_master() && control.print_level("high") && doprint);
   bool oprint = (parall->is_master() && control.print_level("medium") && doprint);
   bool lprint = (parall->is_master() && control.print_level("low") && doprint);

   int ispin = mymolecule.mygrid->ispin;
   int nex[2] = {control.nexcited(0), control.nexcited(1)};
   int neall  = (nex[0] + nex[1]);


   if (neall > 0) 
   {
      coutput << std::endl;
      coutput << " == Virtual Orbital Calculation ==" << std::endl << std::endl;
      coutput << " nex=" << nex[0] << " " << nex[1] << std::endl;
   }
 
}


} // namespace pwdft

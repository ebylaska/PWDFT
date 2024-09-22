
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

/******************************************
 *                                        *
 *            band_cgsd_excited           *
 *                                        *
 ******************************************/
void band_cgsd_excited(Control2 &control, Solid &mysolid, bool doprint, std::ostream &coutput) 
{

   Parallel *parall = mysolid.mygrid->c3db::parall;
   Cneb *mygrid = mysolid.mygrid;
   Ion *myion = mysolid.myion;
 
   bool stalled = false;
   double E[70],total_energy,deltae,deltae_old,deltac;
 
   int it_in   = control.loop(0);
   int it_out  = control.loop(1);
   double tole = control.tolerances(0);
   double tolc = control.tolerances(1);
   double Ep = control.Eprecondition();
   double Sp = control.Sprecondition();
 
   double dt = control.time_step();
   double dte = dt/sqrt(control.fake_mass());
 
   int minimizer   = control.minimizer();
   int lmbfgs_size = control.lmbfgs_size();
 
   bool hprint = (parall->is_master() && control.print_level("high") && doprint);
   bool oprint = (parall->is_master() && control.print_level("medium") && doprint);
   bool lprint = (parall->is_master() && control.print_level("low") && doprint);

   int ispin = mysolid.mygrid->ispin;
   int nex[2] = {control.nexcited(0), control.nexcited(1)};
   int neall  = (nex[0] + nex[1]);


   if (neall > 0) 
   {
      double *vall = mygrid->c_nalloc(ispin);
      // calculating regular virtual orbitals 
      if (lprint) coutput << std::endl;
      if (lprint) coutput << " == Virtual Orbital Calculation ==" << std::endl << std::endl;
      //coutput << " nex=" << nex[0] << " " << nex[1] << std::endl;

      mysolid.epsi_initialize(control.input_e_movecs_filename(),
                              control.input_movecs_initialize(),nex,coutput);

      mysolid.gen_vall();
      mysolid.get_vall(vall);
      mysolid.epsi_minimize(vall,coutput);

      if (lprint) coutput << mysolid.print_virtual();

      //std::cout << "start the exicited minimizer" << std::endl;

    //  std::cout << "check the virtual ortho" << std::endl;
    //  mygrid->g_ortho_excited(mymolecule.psi1,nex, mymolecule.psi1_excited);

    //  std::cout << "calculate the virtual gradient" << std::endl;

    //  std::cout << "task a virtual" << std::endl;
    //  std::cout << "recheck the virtual ortho" << std::endl;

      mysolid.epsi_finalize(control.input_e_movecs_filename(),coutput);

      // read in excited wavefunctions and initialize epsi
     //       call control_ispin_set(psi_ispin())
     //       if (.not.control_check_number_virtuals()) then
     //         call epsi_new()
     //         newpsi = .true.
     //       else
     //         newpsi = .false.

     mygrid->r_dealloc(vall);


   }
 
}

} // namespace pwdft

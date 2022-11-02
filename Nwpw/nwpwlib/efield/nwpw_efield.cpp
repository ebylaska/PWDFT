
#include        <iostream>
#include        <cstring>

#include        "iofmt.hpp"
#include	"nwpw_efield.hpp"


namespace pwdft {


/* Constructor */ 

/*************************************
 *                                   *
 *      nwpw_efield::nwpw_efield     *
 *                                   *
 *************************************/
nwpw_efield::nwpw_efield(Ion *myion0, Pneb *mypneb0, Strfac *mystrfac0, Control2& control, std::ostream& coutput)
{
   mypneb   = mypneb0;
   myion    = myion0;
   mystrfac = mystrfac0;

   ispin  = mypneb->ispin;
   n2ft3d = mypneb->n2ft3d;
   ne    = mypneb->ne;

   double omega = mypneb->lattice->omega();
   double scal1 = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
   dv = omega*scal1;
   bool oprint = ((mypneb->PGrid::parall->is_master()) && (control.print_level("medium")));

   efield_on   = control.Efield_on();

   if (efield_on)  
   {
      efield_type = control.Efield_type();
      efield_vector[0] = control.Efield_vector(0);
      efield_vector[1] = control.Efield_vector(1);
      efield_vector[2] = control.Efield_vector(2);
      efield_center[0] = control.Efield_center(0);
      efield_center[1] = control.Efield_center(1);
      efield_center[2] = control.Efield_center(2);

      v_field =  mypneb->r_alloc();

      double *r_sym_grid = mypneb->r_nalloc(3);
      mypneb->generate_r_sym_grid(r_sym_grid);
      for (auto k=0; k<n2ft3d; ++k)
         v_field[k] +=  efield_vector[0]*(r_sym_grid[3*k]  - efield_center[0])
                     +  efield_vector[1]*(r_sym_grid[3*k+1]- efield_center[1])
                     +  efield_vector[2]*(r_sym_grid[3*k+2]- efield_center[2]);
      mypneb->r_dealloc(r_sym_grid);
   }

   /* write out efield header */
   if (oprint && efield_on)
   {
         coutput << std::endl;
         coutput << " initializing nwpw_efield object" << std::endl;
         coutput << " -------------------------------" << std::endl;
   }

}

/* Object functions */

/***********************************************
 *                                             *
 *        nwpw_efield::shortprint_efield       *
 *                                             *
 ***********************************************/

std::string nwpw_efield::shortprint_efield()
{
   double mu = std::sqrt(efield_vector[0]*efield_vector[0] + efield_vector[1]*efield_vector[1] + efield_vector[2]*efield_vector[2]);
   std::stringstream stream;

   stream << std::endl;
   //stream << "== Electric Field ==" << std::endl << std::endl;
   if (efield_type==0) stream << " periodic Electric field:" << std::endl;
   if (efield_type==1) stream << " APC Electric field:" << std::endl;
   if (efield_type==2) stream << " real space Electric field:" << std::endl;
   stream << "     Electric Field (au) = (" << Ffmt(10,5) << efield_vector[0] << " " 
                                            << Ffmt(10,5) << efield_vector[1] << " " 
                                            << Ffmt(10,5) << efield_vector[2] << " )" << std::endl;
   stream << "             Center (au) = (" << Ffmt(10,5) << efield_center[0] << " " 
                                            << Ffmt(10,5) << efield_center[1] << " " 
                                            << Ffmt(10,5) << efield_center[2] << " )" << std::endl;
   stream << std::endl;

   return stream.str();
}


}

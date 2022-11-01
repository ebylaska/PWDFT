
#include        "iofmt.hpp"
#include	"nwpw_efield.hpp"

namespace pwdft {


/* Constructor */ 

/*************************************
 *                                   *
 *      nwpw_efield::nwpw_efield     *
 *                                   *
 *************************************/
nwpw_efield::nwpw_efield(Ion *myion0, Pneb *mypneb0, Strfac *mystrfac0, Control2&)
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

}

/* Object functions */

/***********************************************
 *                                             *
 *        nwpw_efield::shortprint_efield       *
 *                                             *
 ***********************************************/

std::string nwpw_efield::shortprint_efield()
{
   double mu = std::sqrt(mdipole[0]*mdipole[0] + mdipole[1]*mdipole[1] + mdipole[2]*mdipole[2]);
   std::stringstream stream;

   stream << "== Center of Charge ==" << std::endl << std::endl;
   stream << "spin up    = (" << Ffmt(10,4) << mcdv1[0] << " " 
                              << Ffmt(10,4) << mcdv1[1] << " " 
                              << Ffmt(10,4) << mcdv1[2] << " )" << std::endl;
   stream << "spin down  = (" << Ffmt(10,4) << mcdv2[0] << " " 
                              << Ffmt(10,4) << mcdv2[1] << " " 
                              << Ffmt(10,4) << mcdv2[2] << " )" << std::endl;
   stream << "     total = (" << Ffmt(10,4) << mcdv3[0] << " " 
                              << Ffmt(10,4) << mcdv3[1] << " " 
                              << Ffmt(10,4) << mcdv3[2] << " )" << std::endl;
   stream << "ionic      = (" << Ffmt(10,4) << mqv1[0] << " " 
                              << Ffmt(10,4) << mqv1[1] << " " 
                              << Ffmt(10,4) << mqv1[2] << " )" << std::endl;
   stream << std::endl;
   stream << "== Molecular Dipole wrt Center of Mass ==" << std::endl << std::endl;
   stream << "mu   = (" << Ffmt(10,4) << mdipole[0] << " " 
                        << Ffmt(10,4) << mdipole[1] << " " 
                        << Ffmt(10,4) << mdipole[2] << " ) au" << std::endl;
   stream << "|mu| =  " << Ffmt(10,4) <<  mu           << " au ( " 
                        << Ffmt(10,4) <<  mu*autoDebye << " Debye )"  
                        << std::endl;
   stream << std::endl;

   return stream.str();
}


}

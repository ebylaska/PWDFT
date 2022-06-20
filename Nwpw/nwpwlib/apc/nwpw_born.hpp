#ifndef _NWPW_BORN_HPP_
#define _NWPW_BORN_HPP_

// ****************************************************************************
// *                                                                          *
// *       nwpw_born : used to calculate an extended born model using         *
// *                   APC point charges.                                     *
// *                                                                          *
// * Extended Born solvation model:                                           *
// * G.D. Hawkins, C.R. Cramer,D.G. Truhlar (1995),                           *
// *   Pairwise solute descreening of solute charge from a dielectric medium, *
// *   Chem. Phys. Lett., vol. 246, pages 122-129.                            *
// *                                                                          *
// ****************************************************************************

#include        <iomanip>
#include        <iostream>
#include        <string>

#include        "Control2.hpp"
#include        "Ion.hpp"
#include        "Parallel.hpp"


namespace pwdft {
using namespace pwdft;

class nwpw_born {

   Ion      *myion;
   Parallel *myparall;

public:
   bool born_on,born_relax;

   double *vradii,*bradii;
   double dielec;

   /* constructor */
   nwpw_born(Ion *, Parallel *, Control2&);

   /* destructor */
   ~nwpw_born() {
      if (born_on) {
         delete [] vradii;
         delete [] bradii;
      }
   }

   bool   on()     { return born_on; }
   bool   relax()  { return born_relax; }
   double screen() { return (1.0 - 1.0/dielec); }
   double energy(const double *);
   void   fion(const double *, double *);
   void   dVdq(const double *, double *);

   std::string header_print();
   std::string Qprint(const double *);
   std::string shortprint(const double, const double);

   void writejsonstr(std::string&);



  // friend std::ostream& operator<<(std::ostream& os, const nwpw_born& myborn) {
  //    /* using old style c++ formatting */
  //    std::ios init(NULL);
  //    std::init.copyfmt(os);

   //   os << "     extended Born solvation results    " << std::endl;
   //   os << "     -------------------------------    " << std::endl;
   //   os << "     skipped: no gas phase energy " << << std::endl;
   //   os << "     gas phase energy                 = " << std::scientific << std::setw(19) << std::precision(10) << myborn.egas << std::endl;
   //   os << "     sol phase energy                 = " << std::scientific << std::setw(19) << std::precision(10) << myborn.esol << std::endl;
   //   os << "     (electrostatic) solvation energy = " << std::scientific << std::setw(19) << std::precision(10) << myborn.egas-myborn.esol
   //      << "(" << std::fixed << std::setw(8) << std::precision(3) << (myborn.egas-myborn.esol)*27.2116*23.06 << " kcal/mol)" << std::endl;
   //   os << std::endl << std::endl;

    //  os.copyfmt(init);
    //  return os;
  // }


};

}




#endif

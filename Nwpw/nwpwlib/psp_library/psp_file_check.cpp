#include <iostream>

#include "Control2.hpp"
#include "Ion.hpp"
#include "Parallel.hpp"

#include "psp_library.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *           psp_file_check                *
 *                                         *
 *******************************************/
/*
   Determines whether .psp files are present in the permanent dir
*/
void psp_file_check(Parallel *myparall, Ion *myion, Control2 &control,
                    std::ostream &coutput) {
  double zv;
  bool oprint = ((myparall->is_master()) && (control.print_level("medium")));

  psp_library mypsp_library(control);

  /* write out psp library directories being used */
  if (oprint) {
    coutput << std::endl
            << " psp_library: " << mypsp_library.nwpw_libraryps_dir << std::endl
            << std::endl;
    ;
    for (auto const &x : mypsp_library.libraries)
      coutput << " " << x.first << " library " << x.second << std::endl;
    coutput << std::endl;
  }

  /* Check for psp files - generate them if they do not exist */
  for (auto ia = 0; ia < myion->nkatm; ++ia) {
    if (myparall->is_master()) {
      mypsp_library.psp_check(myion->atom(ia), control, &zv, coutput);
    }
    myparall->Brdcst_Values(0, 0, 1, &zv);
    myion->set_zv_psp(ia, zv);
  }

  /* set the total ion charge in control which in turn sets ispin and ne */
  control.set_total_ion_charge(myion->total_zv());
}

} // namespace pwdft

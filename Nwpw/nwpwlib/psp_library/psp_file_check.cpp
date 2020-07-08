#include        <iostream>
using namespace std;

#include        "Parallel.hpp"
#include        "Control2.hpp"
#include        "Ion.hpp"

#include        "psp_library.hpp"


/*******************************************
 *                                         *
 *           psp_file_check                *
 *                                         *
 *******************************************/
/* 
   Determines whether .psp files are present in the permanent dir 
*/
void psp_file_check(Parallel *myparall, Ion *myion, Control2 &control)
{
   double zv;
   psp_library mypsp_library(control);

   /* write out psp library directories being used */
   if (myparall->is_master())
   {
      std::cout << std::endl << " psp_library: " << mypsp_library.nwpw_libraryps_dir << std::endl << std::endl;;
      for (auto const& x : mypsp_library.libraries)
         std::cout << " " << x.first << " library " << x.second << std::endl;
   }

   /* Check for psp files - generate them if they do not exist */
   for (auto ia=0; ia<myion->nkatm; ++ia)
   {
      if (myparall->is_master())
      {
        mypsp_library.psp_check(myion->atom(ia),control,&zv);
        //std::cout << " zv =" << zv << std::endl;
      }
      myparall->Brdcst_Values(0,0,1,&zv);
      myion->set_zv_psp(ia,zv);
   }

    /* set the total ion charge in control which in turn sets ispin and ne */
    control.set_total_ion_charge(myion->total_zv());
}

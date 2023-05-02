#include <iostream>
#include <sys/stat.h>
#include <unistd.h>

#include "Control2.hpp"
#include "Ion.hpp"
#include "Parallel.hpp"

#include "psp_library.hpp"

// Option for C++17
//#include      <filesystem>
// namespace fs = std::filesystem;

// Option for C++ before C++17

// Option for C++ before C++17
namespace fs0 {
inline bool exists(const std::string &filename) {
   struct stat buffer;
   return (stat(filename.c_str(), &buffer) == 0);
}
} // namespace fs


namespace pwdft {

/*******************************************
 *                                         *
 *           directories_check             *
 *                                         *
 *******************************************/
/*
   Determines whether permanent_dir and scratch_dir are present
*/
bool directories_check(Parallel *myparall, Control2 &control, std::ostream &coutput)
{
   bool has_perm = true;
   bool has_scratch = true;
   int i_perm = 0;
   int i_scratch = 0;
   if (myparall->is_master()) {
      has_perm    = fs0::exists(control.permanent_dir_str); 
      has_scratch = fs0::exists(control.scratch_dir_str); 

      if (!(has_perm && has_scratch)) coutput << std::endl;

      if (!has_perm) 
         coutput << " Could not open a file in permanent directory: " 
                 << control.permanent_dir_str << std::endl;
      if (!has_scratch)
         coutput << " Could not open a file in scratch directory: " 
                 << control.scratch_dir_str << std::endl;

      if (!(has_perm && has_scratch)) 
         coutput << " Both permanent and scratch directories not accessible" << std::endl;
      else if (!(has_perm))
         coutput << " Permanent directory not accessible" << std::endl;
      else if (!(has_scratch)) 
         coutput << " Scratch directory not accessible" << std::endl;

   }
   if (has_perm) i_perm = 1;
   myparall->Brdcst_iValue(0,0,&i_perm);
   has_perm = (i_perm==1);

   if (has_scratch) i_scratch = 1;
   myparall->Brdcst_iValue(0,0,&i_scratch);
   has_scratch = (i_scratch==1);

   return (has_perm && has_scratch);
}

/*******************************************
 *                                         *
 *           psp_file_check                *
 *                                         *
 *******************************************/
/*
   Determines whether .psp files are present in the permanent dir
*/
void psp_file_check(Parallel *myparall, Ion *myion, Control2 &control, std::ostream &coutput) 
{
   double zv;
   bool oprint = ((myparall->is_master()) && (control.print_level("medium")));
 
   psp_library mypsp_library(control);

   /* Check for permanent and scratch directories */
   if (oprint) 
   {
       bool has_perm    = fs0::exists(control.permanent_dir_str); 
       bool has_scratch = fs0::exists(control.scratch_dir_str); 

       if (!(has_perm && has_scratch)) coutput << std::endl;

       if (!has_perm) 
          coutput << " Could not open a file in permanent directory: " 
                  << control.permanent_dir_str << std::endl;
       if (!has_scratch)
          coutput << " Could not open a file in scratch directory: " 
                  << control.scratch_dir_str << std::endl;

       if (!(has_perm && has_scratch)) 
          coutput << " Both permanent and scratch directories not accessible" << std::endl;
       else if (!(has_perm))
          coutput << " Permanent directory not accessible" << std::endl;
       else if (!(has_scratch)) 
          coutput << " Scratch directory not accessible" << std::endl;
   }
 
   /* write out psp library directories being used */
   if (oprint) 
   {
      coutput << std::endl
              << " psp_library: " << mypsp_library.nwpw_libraryps_dir << std::endl
              << std::endl;
     
      for (auto const &x : mypsp_library.libraries)
         coutput << " " << x.first << " library " << x.second << std::endl;
      coutput << std::endl;
   }
 
   /* Check for psp files - generate them if they do not exist */
   for (auto ia=0; ia<myion->nkatm; ++ia) 
   {
      if (myparall->is_master()) {
         mypsp_library.psp_check(myion->atom(ia),control,&zv,coutput);
      }
      myparall->Brdcst_Values(0,0,1,&zv);
      myion->set_zv_psp(ia,zv);
   }
 
   /* set the total ion charge in control which in turn sets ispin and ne */
   control.set_total_ion_charge(myion->total_zv());
}

} // namespace pwdft


#include	<string>
#include	"NwpwLibrarypsConfig.hpp"
#include	"psp_library.hpp"

using namespace std;


/* Constructors */

/*******************************************
 *                                         *
 *      psp_library::psp_library           *
 *                                         *
 *******************************************/
psp_library::psp_library(const std::string  dirname)
{
   /* set  the pseudopotential library directory */
   if (dirname.size()>0) 
   {
       nwpw_libraryps_dir = dirname;
   }
   /* Fetch  the pseudopotential library directory */
   else
   {
       nwpw_libraryps_dir = Nwpw_LIBRARYPS_Default;
       if (const char *libraryps0 = std::getenv("NWPW_LIBRARY"))
          nwpw_libraryps_dir = libraryps0;
       else if (const char *libraryps0 = std::getenv("NWCHEM_NWPW_LIBRARY"))
          nwpw_libraryps_dir = libraryps0;
   }

   default_library = "pspw_default";

}

psp_library::psp_library(Control2& control)
{
   std::string dirname = control.psp_library_dir;

   if (dirname.size()>0) 
   {
       nwpw_libraryps_dir = dirname;
   }
   /* Fetch  the pseudopotential library directory */
   else
   {
       nwpw_libraryps_dir = Nwpw_LIBRARYPS_Default;
       if (const char *libraryps0 = std::getenv("NWPW_LIBRARY"))
          nwpw_libraryps_dir = libraryps0;
       else if (const char *libraryps0 = std::getenv("NWCHEM_NWPW_LIBRARY"))
          nwpw_libraryps_dir = libraryps0;
   }
   default_library = "pspw_default";

   for (auto const& x : control.psp_libraries)
      libraries[x.first] = x.second;
}


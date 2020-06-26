
#include	<cstdlib>
#include	<cstring>
#include	"NwpwLibrarypsConfig.hpp"
#include	"psp_library.hpp"

using namespace std;


/* Constructors */

/*******************************************
 *                                         *
 *      psp_library::psp_library           *
 *                                         *
 *******************************************/
psp_library::psp_library(const char *dirname)
{
   /* set  the pseudopotential library directory */
   if (strlen(dirname)>0) 
   {
       nwpw_libraryps = dirname;
   }
   /* Fetch  the pseudopotential library directory */
   else
   {
       nwpw_libraryps = Nwpw_LIBRARYPS_Default;
       if (const char *libraryps0 = std::getenv("NWPW_LIBRARY"))
          nwpw_libraryps = libraryps0;
       else if (const char *libraryps0 = std::getenv("NWCHEM_NWPW_LIBRARY"))
          nwpw_libraryps = libraryps0;
   }

}


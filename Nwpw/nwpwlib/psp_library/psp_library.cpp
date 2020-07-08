
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;

#include	<string>
#include	"NwpwLibrarypsConfig.hpp"
#include	"psp_library.hpp"

using namespace std;


static int convert_psp_type(char *test)
{
   int psp_type = 0;
   if (test[0]=='0') psp_type = 0;
   if (test[0]=='1') psp_type = 1;
   if (test[0]=='2') psp_type = 2;
   if (test[0]=='3') psp_type = 3;
   if (test[0]=='4') psp_type = 4;
   if (test[0]=='5') psp_type = 5;
   if (test[0]=='6') psp_type = 6;
   if (test[0]=='7') psp_type = 7;
   if (test[0]=='8') psp_type = 8;
   if (test[0]=='9') psp_type = 9;

   return psp_type;
}

/*******************************************
 *                                         *
 *            psp_read_header              *
 *                                         *
 *******************************************/

/* This function returns the header data of a psp file.  

   Entry - myparall
           fname - name of vpp file
   Exit -  header data: comment, psp_type, version, nfft, unita, atom, amass, zv

   Returns true if fname exists otherwise false.
*/

#define FMT1    "%lf"

static bool psp_read_header(char *fname, double *zv)
{
   int i,ifound,ihasae;
   char atom[2];
   FILE *fp;

   ifound = cfileexists(fname);


   if (ifound>0)
   {
      fp = std::fopen(fname,"r");
      std::fscanf(fp,"%s",atom);
      ihasae = convert_psp_type(atom);
      if (ihasae>0) std::fscanf(fp,"%s",atom);
      
      std::fscanf(fp,FMT1,zv);
      std::fclose(fp);
   }
   return (ifound>0);
}



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


/*******************************************
 *                                         *
 *      psp_library::psp_check             *
 *                                         *
 *******************************************/

void psp_library::psp_check(const char *atom, Control2& control, double *zv)
{

   char fname[256];

   strcpy(fname,atom);
   strcat(fname,".psp");
   control.add_permanent_dir(fname);
   if (!psp_read_header(fname,zv))
   {
      std::cout << " -- need to generate psp =" << fname << std::endl;

      *zv = 0.0;
   }
}

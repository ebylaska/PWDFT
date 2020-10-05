
#include	<string>
#include	<cstring>
#include	<iostream>
#include	<cstdio>
#include	<cstdlib>
#include	<cmath>

#include	"parsestring.hpp"

extern "C" {
#include        "pseudopotential.h"
}


using namespace std;


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
   nwpw_permanent_dir = "";

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
   default_library    = "pspw_default";
   nwpw_permanent_dir = string(control.permanent_dir());

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
   strcat(fname,".psp"); control.add_permanent_dir(fname);

   if (!cfileexists(fname))
      this->psp_generator_auto(atom,control);

   if (!psp_read_header(fname,zv))
   {
      std::cout << " -- need to generate psp =" << fname << std::endl;

      *zv = 0.0;
   }
}


/*******************************************
 *                                         *
 *     psp_library::psp_generator_auto     *
 *                                         *
 *******************************************/
void psp_library::psp_generator_auto(const char *atom, Control2& control)
{
   char fname[256],tname[256];
   char sdir_name[256], dir_name[256];
   char psp_in[10],psp_out[10];
   char aatom[10];

   int ptype = this->psp_type(atom);
   int lmax = this->psp_lmax(atom);
   int locp = this->psp_locp(atom);
   double rlocal = this->psp_rlocal(atom);
   //std::cout << " psptype = " << ptype << std::endl;
   //std::cout << " lmax = " << lmax << std::endl;
   //std::cout << " locp = " << locp << std::endl;
   //std::cout << " rlocal = " << rlocal << std::endl;

   strcpy(aatom,atom);
   strcpy(fname,atom);
   strcat(fname,".psp");
   strcpy(psp_out,fname); 
   control.add_permanent_dir(fname);
   std::string mypsp_fname(fname);
   std::string mypsp_libname = this->psp_libname(atom);

   strcpy(tname,"junk");
   strcat(tname,".inp"); 
   strcpy(psp_in,tname); 
   control.add_permanent_dir(tname);
   std::string mypsp_junkname(tname);

   std::string pspinput = mystring_readfile(mypsp_libname);

   mystring_writefile(mypsp_junkname,pspinput);

   std::cout << " generating 1d pseudopotential file: " << mypsp_fname << std::endl;

   strcpy(sdir_name,control.scratch_dir());
   strcpy(dir_name,control.permanent_dir());
   int n9 = strlen(sdir_name);
   int n0 = strlen(dir_name);
   int n1 = strlen(psp_in);
   int n2 = strlen(psp_out);
   int n3 = strlen(aatom);
   int efg00 = 0;

   int   debug = 0;
   int   print = 0;
   debug = control.print_level("debug");
   print = control.print_level("high");

   if ((ptype==0)||(ptype==7))
   {
      //std::cout << " sdir_name = " << sdir_name << " " << strlen(sdir_name) << std::endl;
      //std::cout << " dir_name = " << dir_name << " " << strlen(dir_name) << std::endl;
      //std::cout << " psp_in = " << psp_in << " " << strlen(psp_in) << std::endl;
      //std::cout << " psp_out = " << psp_out << " " << strlen(psp_out) << std::endl;
      pspsolve(&print,&debug,&lmax,&locp,&rlocal,&efg00,
               sdir_name,&n9,
               dir_name,&n0,
               psp_in,&n1,
               psp_out,&n2);

   }
   else if (ptype==1)
   {
      std::cout << " -- psp_type=1 not tested" << std::endl;
      hgh_parse(&debug,&lmax,&locp,&rlocal,
                sdir_name,&n9,
                dir_name,&n0,
                psp_in,&n1,
                psp_out,&n2);

   }
   else if (ptype==2)
   {
      std::cout << " -- psp_type=2 not tested" << std::endl;
      cpi_parse(&debug,&lmax,&locp,&rlocal,
                sdir_name,&n9,
                dir_name,&n0,
                psp_in,&n1,
                psp_out,&n2,
                aatom,&n3);
   }
   else if (ptype==3)
   {
      std::cout << " -- psp_type=3 not tested" << std::endl;
      teter_parse(&debug,&lmax,&locp,&rlocal,
                  sdir_name,&n9,
                  dir_name,&n0,
                  psp_in,&n1,
                  psp_out,&n2,
                  aatom,&n3);
   }
   else if (ptype==4)
   {
      std::cout << " -- psp_type=4 not tested" << std::endl;
      paw_atom_driver(&debug,&lmax,&locp,&rlocal,
                      sdir_name,&n9,
                      dir_name,&n0,
                      psp_in,&n1,
                      psp_out,&n2,
                      aatom,&n3);
   }
   else if (ptype==5)
   {
      std::cout << " -- psp_type=5 not tested" << std::endl;
      qmmm_parse(&debug,&lmax,&locp,&rlocal,
                 sdir_name,&n9,
                 dir_name,&n0,
                 psp_in,&n1,
                 psp_out,&n2,
                 aatom,&n3);
   }
   else if (ptype==6)
   {
      std::cout << " -- psp_type=6 not tested" << std::endl;
      carter_parse(&debug,&lmax,&locp,&rlocal,
                   sdir_name,&n9,
                   dir_name,&n0,
                   psp_in,&n1,
                   psp_out,&n2,
                   aatom,&n3);
   }
   else
   {
      std::cout << " -- unknown psp_type, psp_type=" << ptype << std::endl;
   }

   /* delete the junk.inp file */
  std::remove(tname);

}



/*******************************************
 *                                         *
 *        psp_library::psp_libname         *
 *                                         *
 *******************************************/
string psp_library::psp_libname(const char *atom)
{

   /* filename of psp file for atom */
   std::string aname(atom);
   std::string libname = nwpw_libraryps_dir + "/" ;
   if (libraries.count(aname)>0)
      libname += libraries[aname] + "/" + aname;
   else
      libname += default_library + "/" + aname;

   return libname;
}

/*******************************************
 *                                         *
 *        psp_library::psp_type            *
 *                                         *
 *******************************************/

int psp_library::psp_type(const char *atom)
{
   /* filename of psp file for atom */
   /* load the file into string */
   std::string aa = mystring_readfile(this->psp_libname(atom));

   int ptype = 0;
   if (mystring_contains(aa,"[psp_type]"))
      ptype = std::stoi(mystring_trim(mystring_split(mystring_split(aa,"[psp_type]")[1],"\n")[0]));

   return ptype;
}

/*******************************************
 *                                         *
 *        psp_library::psp_lmax            *
 *                                         *
 *******************************************/
int psp_library::psp_lmax(const char *atom)
{
   /* filename of psp file for atom */
   /* load the file into string */
   std::string aa = mystring_readfile(this->psp_libname(atom));

   int lmax = -1;
   if (mystring_contains(aa,"[lmax]"))
      lmax = std::stoi(mystring_trim(mystring_split(mystring_split(aa,"[lmax]")[1],"\n")[0]));

   return lmax;
}

/*******************************************
 *                                         *
 *        psp_library::psp_locp            *
 *                                         *
 *******************************************/
int psp_library::psp_locp(const char *atom)
{
   /* filename of psp file for atom */
   /* load the file into string */
   std::string aa = mystring_readfile(this->psp_libname(atom));

   int locp = -1;
   if (mystring_contains(aa,"[locp]"))
      locp = std::stoi(mystring_trim(mystring_split(mystring_split(aa,"[locp]")[1],"\n")[0]));

   return locp;
}

/*******************************************
 *                                         *
 *        psp_library::psp_rlocal          *
 *                                         *
 *******************************************/
double psp_library::psp_rlocal(const char *atom)
{
   /* filename of psp file for atom */
   /* load the file into string */
   std::string aa = mystring_readfile(this->psp_libname(atom));

   double rlocal = 1.0;
   if (mystring_contains(aa,"[rlocal]"))
      rlocal = std::stod(mystring_trim(mystring_split(mystring_split(aa,"[rlocal]")[1],"\n")[0]));

   return rlocal;
}





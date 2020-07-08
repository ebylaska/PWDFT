#ifndef _PSP_LIBRARY_HPP_
#define _PSP_LIBRARY_HPP_
/* psp_library.h -
   author - Eric Bylaska

*/

//extern "C" {
//#include        "compressed_io.h"
//}
#include        "compressed_io.hpp"



#include	<map>
#include	<string>
#include	"Control2.hpp"

class psp_library {

   public:
   std::string nwpw_libraryps_dir;
   std::string default_library;
   std::map<std::string, std::string> libraries;

   /* Constructors */
   psp_library(const std::string);
   psp_library(Control2&);

   /* other */
   void psp_check(const char *, Control2&, double *);

};


#endif

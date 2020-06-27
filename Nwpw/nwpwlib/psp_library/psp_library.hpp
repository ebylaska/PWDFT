#ifndef _PSP_LIBRARY_HPP_
#define _PSP_LIBRARY_HPP_
/* psp_library.h -
   author - Eric Bylaska

*/
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

};


#endif

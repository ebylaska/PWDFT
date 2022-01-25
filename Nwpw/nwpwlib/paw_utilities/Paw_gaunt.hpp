#ifndef _PAW_GAUNT_HPP_
#define _PAW_GAUNT_HPP_

namespace pwdft {
using namespace pwdft;


class Paw_gaunt {

public:
   bool gaunt_iscmpx;
   int  gaunt_lmax;
   int  gaunt_sizel2;
   double *gaunt_coeff;

   /* Constructors */
   Paw_gaunt(const bool, const int);

   /* destructor */
   ~Paw_gaunt() {
      delete [] gaunt_coeff;
   }

   /* Gaunt routines */
   double gaunt(const int, const int, const int, const int, const int, const int);
   double gaunt2(const int, const int, const int, const int, const int, const int);
   double gaunt3(const int, const int, const int, const int, const int, const int);
};

}

#endif

#ifndef	_STRFAC_H_
#define _STRFAC_H_



#include        "PGrid.hpp"
#include	"Ion.hpp"

class	Strfac {
   int    *i_indx[2],*j_indx[2],*k_indx[2];
   double *wx1,*wy1,*wz1;
   double unita[9],unitg[9];

   PGrid  *mygrid;
   Ion	  *myion;
public:

   /* Constructors */
   Strfac(Ion *, PGrid *);

   /* destructor */
   ~Strfac() {
            delete [] i_indx[0];
            delete [] j_indx[0];
            delete [] k_indx[0];
            delete [] i_indx[1];
            delete [] j_indx[1];
            delete [] k_indx[1];
            delete [] wx1;
            delete [] wy1;
            delete [] wz1;
         }

    void phafac();
    void strfac_pack(const int, const int, double *);

};

#endif

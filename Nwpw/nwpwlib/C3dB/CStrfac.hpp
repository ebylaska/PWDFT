#ifndef _CSTRFAC_H_
#define _CSTRFAC_H_

#pragma once

#include "Ion.hpp"
#include "CGrid.hpp"

namespace pwdft {

class CStrfac {
  //int *i_indx[2], *j_indx[2], *k_indx[2];
  int maxsize;
  int **i_indx, **j_indx, **k_indx;
  double *wx1, *wy1, *wz1;
  double unita[9], unitg[9];

  CGrid *mygrid;
  Ion *myion;

public:
   /* Constructors */
   CStrfac(Ion *, CGrid *);
 
   /* destructor */
   ~CStrfac() {
      for (auto nb=0; nb<maxsize; ++nb)
      {
         delete[] i_indx[nb];
         delete[] j_indx[nb];
         delete[] k_indx[nb];
      }
      delete[] i_indx;
      delete[] j_indx;
      delete[] k_indx;
 
      delete[] wx1;
      delete[] wy1;
      delete[] wz1;
   }
 
   void phafac();
   void strfac_pack(const int, const int, double *);
};
} // namespace pwdft

#endif

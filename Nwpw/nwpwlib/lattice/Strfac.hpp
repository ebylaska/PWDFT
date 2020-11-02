#ifndef	_STRFAC_H_
#define _STRFAC_H_

#pragma once

#include        "PGrid.hpp"
#include	"Ion.hpp"
// NOTE: the SYCL segment is commented out for the lack
// of memory errors
// #ifdef NWPW_SYCL
// #include "gdevice.hpp"
// #endif

class	Strfac {
   int    *i_indx[2],*j_indx[2],*k_indx[2];
   double *wx1,*wy1,*wz1;
// #ifdef NWPW_SYCL
//    int *i1_indx_sycl, *j1_indx_sycl, *k1_indx_sycl;
//    double *wx1_sycl, *wy1_sycl, *wz1_sycl;
// #endif
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
// #ifdef NWPW_SYCL
// 	    cl::sycl::free(wx1_sycl, *get_syclQue());
// 	    cl::sycl::free(wy1_sycl, *get_syclQue());
// 	    cl::sycl::free(wz1_sycl, *get_syclQue());

// 	    cl::sycl::free(i1_indx_sycl, *get_syclQue());
// 	    cl::sycl::free(j1_indx_sycl, *get_syclQue());
// 	    cl::sycl::free(k1_indx_sycl, *get_syclQue());
// #endif
         }

    void phafac();
    void strfac_pack(const int, const int, double *);
// #ifdef NWPW_SYCL
//     void strfac_pack_sycl(const int, const int, double *);
// #endif
};

#endif

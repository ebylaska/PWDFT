#ifndef _nwpw_lmbfgs_HPP_
#define _nwpw_lmbfgs_HPP_

#pragma once

/* nwpw_lmbfgs.hpp
   Author - Eric Bylaska
        this class is used to keep perform lmfgs steps


*/

#include	<cmath>
#include	"blas.h"

namespace pwdft {


class nwpw_lmbfgs {

   int one  = 1;
   double rone  = 1.0;
   double mrone = -1.0;

   int    m,max_m,nsize;
   int    *indx;
   double *rho;
   double *sylist;

public: 

   /* constructor */
   nwpw_lmbfgs(int nsize0, int max_m0, double *x0, double *g0) {
      m     = 0;
      max_m = max_m0;
      nsize = nsize0;

      indx   = new int [max_m];
      rho    = new double [max_m];
      sylist = new double [2*(max_m)*nsize];

      for (auto k=0; k<max_m; ++k) indx[k] = k;

      //DCOPY_PWDFT(nsize,x0,one,&sylist[2*m*nsize],    one);
      //DCOPY_PWDFT(nsize,g0,one,&sylist[(2*m+1)*nsize],one);
      std::memcpy(sylist+(2*m)  *nsize,x0,nsize*sizeof(double));
      std::memcpy(sylist+(2*m+1)*nsize,g0,nsize*sizeof(double));
   }

   /* destructor */
   ~nwpw_lmbfgs() {
       delete [] sylist;
       delete [] rho;
       delete [] indx;
   }

   void lmbfgs(double *x, double *g, double *q) {
      //DCOPY_PWDFT(nsize,g,one,q,one);
      //DCOPY_PWDFT(nsize,x,one,&sylist[(2*(indx[m+1]))  *nsize],one);
      //DCOPY_PWDFT(nsize,g,one,&sylist[(2*(indx[m+1])+1)*nsize],one);
      std::memcpy(q,g,nsize*sizeof(double));
      std::memcpy(sylist+(2*(indx[m+1]))  *nsize,x,nsize*sizeof(double));
      std::memcpy(sylist+(2*(indx[m+1])+1)*nsize,g,nsize*sizeof(double));

      DAXPY_PWDFT(nsize,mrone,x,one,&sylist[(2*indx[m])  *nsize],one);
      DAXPY_PWDFT(nsize,mrone,g,one,&sylist[(2*indx[m]+1)*nsize],one);

      double sum  = DDOT_PWDFT(nsize,&sylist[(2*indx[m]+1)*nsize],one,&sylist[(2*indx[m])*nsize],one);
      if (std::fabs(sum)>1.0e-11) {
         rho[indx[m]] = 1.0/sum;

         double alpha[m+1];
         for (auto k=m-1; k>=0; --k) {
            alpha[k]     = rho[indx[k]]*DDOT_PWDFT(nsize,&sylist[(2*indx[k])*nsize],one,q,one);
            double tscal = -alpha[k];
            DAXPY_PWDFT(nsize,tscal,&sylist[(2*indx[k]+1)*nsize],one,q,one);
         }

         /* add preconditioner here */
         if (m>200) {
             double sumsy = DDOT_PWDFT(nsize,&sylist[(2*indx[m-1])*nsize],  one,&sylist[(2*indx[m-1]+1)*nsize],one);
             double sumyy = DDOT_PWDFT(nsize,&sylist[(2*indx[m-1]+1)*nsize],one,&sylist[(2*indx[m-1]+1)*nsize],one);
             if (sumyy>1.0e-10) {
                double gmma = sumsy/sumyy;
                if (gmma>2.0) gmma = 2.0;
                DSCAL_PWDFT(nsize,gmma,q,one);
             }
         }

         for (auto k=0; k<m; ++k) {
            double beta  = rho[indx[k]]*DDOT_PWDFT(nsize,&sylist[(2*indx[k]+1)*nsize],one,q,one);
            double tscal = -(beta-alpha[k]);
            DAXPY_PWDFT(nsize,tscal,&sylist[(2*indx[k])*nsize],one,q,one);
         }

         if (m<(max_m-2)) {
            ++m;
         } else {
            int itmp = indx[0];
            for (auto k=0; k<(max_m-1); ++k) 
               indx[k]=indx[k+1];
            indx[max_m-1] = itmp;
         }
      }

   }

};

}

#endif

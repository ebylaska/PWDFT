#ifndef _nwpw_ckerker_HPP_
#define _nwpw_ckerker_HPP_

#pragma once

/* nwpw_ckerker.hpp
   Author - Eric Bylaska
        this class is used to keep perform kerker updates


*/

#include "blas.h"
#include <cmath>
#include "CGrid.hpp"

namespace pwdft {

class nwpw_ckerker {

   CGrid *mygrid;
   bool dokerker = false;
   bool isband = false;;
   int npack0, n2ft3d;
   double *tg,*tmpv;
   double scal1,g0;

public:
   /* constructors */
   /*******************************************
    *                                         *
    *        nwpw_ckerker::nwpw_ckerker       *
    *                                         *
    *******************************************/
   nwpw_ckerker(CGrid *mygrid0, const double ing0) 
   {
      mygrid = mygrid0;
      g0 = ing0;
      isband = false;
      dokerker = (g0>0.0);
      n2ft3d     = mygrid->n2ft3d;
      npack0     = mygrid->npack(0);
      if (dokerker)
      {
         double gg0 = g0*g0;
         double *Gpackx = mygrid->Gpackxyz(0,0);
         double *Gpacky = mygrid->Gpackxyz(0,1);
         double *Gpackz = mygrid->Gpackxyz(0,2);
         tg   = new double[npack0];
         tmpv = new double[n2ft3d];
         for (auto k=0; k<npack0; ++k)
         {
            double gx = Gpackx[k];
            double gy = Gpacky[k];
            double gz = Gpackz[k];
            double gg = gx*gx + gy*gy + gz*gz;
            tg[k] = gg/(gg+gg0);
         }
      }
      scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
   }



   /* destructor */
   /*******************************************
    *                                         *
    *        nwpw_ckerker::~nwpw_ckerker      *
    *                                         *
    *******************************************/
   ~nwpw_ckerker() 
   {
      if (dokerker)
      {
         delete[] tg;
         delete[] tmpv;
      }
   }

   /*******************************************
    *                                         *
    *           nwpw_ckerker::kerker_G        *
    *                                         *
    *******************************************/
   void kerker_G(double *v) 
   {
      if (dokerker)
      {
         mygrid->rc_SMul(scal1,v,tmpv);
         mygrid->rc_pfft3f(0,tmpv);
         mygrid->c_pack(0,tmpv);
       
         int k1=0;
         int k2=0;
         for (auto k=0; k<npack0; ++k)
         {
            tmpv[k1] = tmpv[k1]*tg[k];
            tmpv[k2] = tmpv[k2]*tg[k];
            k1 += 2;
            k2 += 2;
         }
         mygrid->c_unpack(0,tmpv);
         mygrid->cr_pfft3b(0,tmpv);
         mygrid->cr_copy(tmpv,v);
      }
   }
};

} // namespace pwdft

#endif
